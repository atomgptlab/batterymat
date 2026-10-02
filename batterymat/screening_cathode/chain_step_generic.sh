#!/bin/bash
#SBATCH --job-name=bmchain
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=chain_%x_%j.out
#SBATCH --error=chain_%x_%j.err

# =============================================================================
# Cluster-generic copy of chain_step.sh (self-resubmitting delithiation chain).
# EDIT THE BLOCK BELOW for your cluster, and the #SBATCH header above
# (partition, --account=..., ntasks, memory, walltime). #SBATCH lines cannot
# use shell variables, so they must be edited by hand.
#
# Usage (from $WORKDIR, which must contain dft_prep.py and dft_inputs/):
#   sbatch --job-name=bm_JVASP-2017-LCO-PBE  chain_step_generic.sh JVASP-2017-LCO-PBE
#   sbatch --job-name=bm_JVASP-144791-NMC-vdW chain_step_generic.sh JVASP-144791-NMC-vdW
#   sbatch --job-name=bm_JVASP-79809-NCO      chain_step_generic.sh JVASP-79809-NCO
#   sbatch --job-name=bm_JVASP-11340-MMO      chain_step_generic.sh JVASP-11340-MMO
# The argument is the DIRECTORY NAME under dft_inputs/ (not the bare JID:
# JVASP-2017 is ambiguous now that JVASP-2017-LCO and JVASP-2017-LCO-PBE coexist).
#
# Each invocation: runs VASP on the newest step, copies outputs to results/,
# records the energy (auto-read from results/OUTCAR), generates the next step
# (ALIGNN vacancy ranking, needs the alignn conda env), builds its POTCAR, and
# resubmits itself. When no further step can be generated it computes the
# voltage curve and stops.
# =============================================================================

# ----------------------------- EDIT ME ---------------------------------------
WORKDIR="${BMAT_WORKDIR:-/data/$USER/revision}"   # contains dft_prep.py + dft_inputs/
CONDA_SH="${BMAT_CONDA_SH:-/data/$USER/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${BMAT_CONDA_ENV:-alignn}"                 # env with numpy/jarvis-tools/alignn/torch
VASP_MODULE="${BMAT_VASP_MODULE:-}"                   # e.g. "vasp/6.4.3"; leave empty if none
VASP_BIN="${BMAT_VASP_BIN:-}"                         # full path to vasp_std; empty = auto-detect
POTCAR_DIR="${BMAT_POTCAR_DIR:-${VASP_PP_PATH:-}}"    # dir containing Li_sv/POTCAR, Na_pv/POTCAR, ...
VDW_KERNEL="${BMAT_VDW_KERNEL:-}"                     # path to vdw_kernel.bindat (optB88-vdW runs);
                                                      # empty = VASP builds it on first use (slow but fine)
E_ION_METAL="${BMAT_E_ION_METAL:-}"                   # Na/Mg chains: metal reference (eV/atom) once
                                                      # ref-Na / ref-Mg have run; empty = skip voltage
# -----------------------------------------------------------------------------

set -e
set -o pipefail

JID="$1"
[ -z "$JID" ] && { echo "usage: sbatch chain_step_generic.sh <dft_inputs dir name>"; exit 1; }
# Scavenger jobs can be preempted and requeued; if another segment of this chain is
# already running (requeue raced a resubmission), stop so two jobs never share a step.
OTHERS=$( { squeue -h -u "$USER" -n "bm_$JID" -t RUNNING -o %i 2>/dev/null || true; } | { grep -vx "${SLURM_JOB_ID:-x}" || true; } | wc -l )
if [ "$OTHERS" -gt 0 ]; then echo "== $JID: another running segment found, exiting"; exit 0; fi
cd "$WORKDIR"

source "$CONDA_SH"
conda activate "$CONDA_ENV"
[ -d "${ALIGNN_REPO:-$HOME/alignn/repo}/alignn" ] && export PYTHONPATH="${ALIGNN_REPO:-$HOME/alignn/repo}${PYTHONPATH:+:$PYTHONPATH}"   # skipjack only
export OMP_NUM_THREADS=1      # hybrid MPI+OpenMP VASP builds otherwise oversubscribe the node
ulimit -s unlimited

[ -n "$VASP_MODULE" ] && module load "$VASP_MODULE"
if [ -z "$VASP_BIN" ]; then
    command -v vasp_std >/dev/null 2>&1 && VASP_BIN=$(command -v vasp_std)
fi
[ -z "$VASP_BIN" ] && { echo "FATAL: vasp_std not found. Set VASP_BIN at the top of this script."; exit 1; }
[ -z "$POTCAR_DIR" ] || [ ! -d "$POTCAR_DIR" ] && { echo "FATAL: POTCAR_DIR not set/found. Edit the top of this script."; exit 1; }
echo "Using VASP:    $VASP_BIN (ranks scaled per step, max $SLURM_NTASKS)"
echo "Using POTCARs: $POTCAR_DIR"

SUP_DIR=$(ls -d dft_inputs/"$JID"/supercell_* | head -1)
STEP_DIR=$(ls -d "$SUP_DIR"/step_* | sort -t_ -k2 -n | tail -1)
STEP_NUM=$(basename "$STEP_DIR" | sed -E 's/step_0*([0-9]+)_.*/\1/')

build_potcar () {
    # POTCAR_spec lists one PAW label per line (e.g. Li_sv, Na_pv, Mg_pv, Co, O)
    local d="$1"
    [ -s "$d/POTCAR" ] && return 0
    : > "$d/POTCAR"
    while read -r el; do
        [ -z "$el" ] && continue
        cat "$POTCAR_DIR/$el/POTCAR" >> "$d/POTCAR"
    done < "$d/POTCAR_spec"
    # optB88-vdW runs need vdw_kernel.bindat in the run directory (or VASP generates it)
    if grep -qi "LUSE_VDW" "$d/INCAR" && [ -n "$VDW_KERNEL" ] && [ ! -f "$d/vdw_kernel.bindat" ]; then
        cp "$VDW_KERNEL" "$d/vdw_kernel.bindat"
    fi
}

if [ ! -f "$STEP_DIR/results/OUTCAR" ]; then
    # Scale MPI ranks to the structure and keep NP a multiple of KPAR (INCAR).
    NIONS=$(awk 'NR==7{for(i=1;i<=NF;i++)s+=$i; print s}' "$STEP_DIR/POSCAR")
    KPAR=$(awk -F= 'toupper($0) ~ /^ *KPAR/ {print $2+0}' "$STEP_DIR/INCAR" | head -1)
    KPAR=${KPAR:-1}; [ "$KPAR" -lt 1 ] && KPAR=1
    NP=$SLURM_NTASKS
    [ -n "$NIONS" ] && [ "$NIONS" -lt "$NP" ] && NP=$NIONS
    NP=$(( NP / KPAR * KPAR ))
    [ "$NP" -lt "$KPAR" ] && NP=$KPAR
    # NCORE must divide the ranks per k-point group, or VASP aborts in M_divide on small cells.
    NC=$(awk -F= 'toupper($0) ~ /^ *NCORE/ {print $2+0}' "$STEP_DIR/INCAR" | head -1)
    if [ -n "$NC" ] && [ "$NC" -gt 1 ]; then
        PER=$(( NP / KPAR ))
        while [ $(( PER % NC )) -ne 0 ]; do NC=$(( NC - 1 )); done
        sed -i "s/^ *NCORE *=.*/NCORE = $NC/" "$STEP_DIR/INCAR"
    fi
    echo "== $JID step $STEP_NUM: running VASP with $NP ranks ($NIONS ions, KPAR=$KPAR) in $STEP_DIR ($(date))"
    build_potcar "$STEP_DIR"
    ( cd "$STEP_DIR" && mpirun -np "$NP" "$VASP_BIN" )
    mkdir -p "$STEP_DIR/results"
    cp "$STEP_DIR"/OUTCAR "$STEP_DIR"/CONTCAR "$STEP_DIR"/OSZICAR "$STEP_DIR/results/"
    rm -f "$STEP_DIR"/WAVECAR "$STEP_DIR"/CHG   # never reused by the next step; saves scratch space
    grep -q "reached required accuracy" "$STEP_DIR/results/OUTCAR" || echo "NOT_CONVERGED: $STEP_DIR hit NSW without reaching EDIFFG; energy recorded anyway, inspect before use"
else
    echo "== $JID step $STEP_NUM: results/OUTCAR already present, skipping VASP"
fi

VOLT_ARGS=""
[ -n "$E_ION_METAL" ] && VOLT_ARGS="--e-ion-metal $E_ION_METAL"

echo "== $JID: recording step $STEP_NUM"
python dft_prep.py record "$JID" "$STEP_NUM"

echo "== $JID: generating next step"
if python dft_prep.py next "$JID"; then
    NEW_STEP=$(ls -d "$SUP_DIR"/step_* | sort -t_ -k2 -n | tail -1)
    if [ "$NEW_STEP" = "$STEP_DIR" ]; then
        echo "== $JID: next reported success but produced no new step, stopping to avoid a loop"
        python dft_prep.py voltage "$JID" $VOLT_ARGS || true
        exit 0
    fi
    build_potcar "$NEW_STEP"
    echo "== $JID: resubmitting chain for $(basename "$NEW_STEP")"
    sbatch --job-name="bm_$JID" "$0" "$JID"
else
    echo "== $JID: no further step (delithiation complete), computing voltage curve"
    python dft_prep.py voltage "$JID" $VOLT_ARGS || true
fi
echo "== $JID: chain segment done ($(date))"
