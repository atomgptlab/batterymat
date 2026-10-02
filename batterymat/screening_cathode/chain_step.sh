#!/bin/bash
#SBATCH --job-name=bmchain
#SBATCH --account=kchoudh2
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=chain_%x_%j.out
#SBATCH --error=chain_%x_%j.err

# Self-resubmitting delithiation chain for one BatteryMat prospective candidate.
# All #SBATCH directives are above this line on purpose: set -e before them
# silently disables them (same trap noted in job_stage4_relax.sh).
#
# Usage:
#   sbatch --job-name=bm_JVASP-80802 chain_step.sh JVASP-80802     # test on smallest first
#   for j in $(ls dft_inputs); do sbatch --job-name="bm_$j" chain_step.sh "$j"; done
#
# Each invocation: runs VASP on the newest step, copies outputs to results/,
# records the energy (auto-read from results/OUTCAR), generates the next step
# (ALIGNN Li-vacancy ranking), builds its POTCAR, and resubmits itself.
# When no further step can be generated it computes the voltage curve and stops.

set -e
set -o pipefail

JID="$1"
WORKDIR="/data/$USER/dft/neurips"
cd "$WORKDIR"

# --- environment (atomgptlab conda, alignn env has numpy/jarvis/alignn) ---
source /data/$USER/miniforge3/etc/profile.d/conda.sh
conda activate alignn
# vasp_std here is an MPI+OpenMP build: without this it spawns 128 threads per
# rank (2,048 threads on one node) and a 4-atom relax crawls
export OMP_NUM_THREADS=1
# VASP segfaults at startup without a large stack (seen on this cluster)
ulimit -s unlimited

# --- locate VASP (edit VASP_CMD directly if auto-detection picks wrong) ---
if command -v vasp_std >/dev/null 2>&1; then
    VASP_BIN=$(command -v vasp_std)
elif module load vasp 2>/dev/null && command -v vasp_std >/dev/null 2>&1; then
    VASP_BIN=$(command -v vasp_std)
else
    VASP_BIN=$(ls /opt/vasp*/bin/vasp_std /data/apps/vasp*/bin/vasp_std /usr/local/vasp*/bin/vasp_std 2>/dev/null | head -1 || true)
fi
[ -z "$VASP_BIN" ] && { echo "FATAL: vasp_std not found. Set VASP_BIN in this script."; exit 1; }

# --- locate POTCAR library (edit POTCAR_DIR directly if auto-detection fails) ---
POTCAR_DIR="${VASP_PP_PATH:-}"
if [ -z "$POTCAR_DIR" ] || [ ! -d "$POTCAR_DIR" ]; then
    POTCAR_DIR=$(ls -d /data/*/potcar*/potpaw_PBE* /data/*/vasp*/potpaw_PBE* /opt/vasp*/potpaw_PBE* ~/potpaw_PBE* 2>/dev/null | head -1 || true)
fi
[ -z "$POTCAR_DIR" ] && { echo "FATAL: POTCAR library not found. Set POTCAR_DIR in this script."; exit 1; }
echo "Using VASP:    $VASP_BIN (ranks scaled per step, max $SLURM_NTASKS)"
echo "Using POTCARs: $POTCAR_DIR"

SUP_DIR=$(ls -d dft_inputs/"$JID"/supercell_* | head -1)
STEP_DIR=$(ls -d "$SUP_DIR"/step_* | sort -t_ -k2 -n | tail -1)
STEP_NUM=$(basename "$STEP_DIR" | sed -E 's/step_0*([0-9]+)_.*/\1/')

build_potcar () {
    # POTCAR_spec lists one PAW label per line (e.g. Li_sv, Fe_pv, O)
    local d="$1"
    [ -s "$d/POTCAR" ] && return 0
    : > "$d/POTCAR"
    while read -r el; do
        [ -z "$el" ] && continue
        cat "$POTCAR_DIR/$el/POTCAR" >> "$d/POTCAR"
    done < "$d/POTCAR_spec"
}

if [ ! -f "$STEP_DIR/results/OUTCAR" ]; then
    # Scale MPI ranks to the structure: more ranks than ions over-decomposes
    # tiny cells (16 ranks on 4 ions segfaulted at startup on this cluster).
    # Rank count must also be a multiple of KPAR from the INCAR: 3 ranks with
    # KPAR=2 aborts in M_divide (seen on step_01_Li0 of the test chain).
    NIONS=$(awk 'NR==7{for(i=1;i<=NF;i++)s+=$i; print s}' "$STEP_DIR/POSCAR")
    KPAR=$(awk -F= 'toupper($0) ~ /^ *KPAR/ {print $2+0}' "$STEP_DIR/INCAR" | head -1)
    KPAR=${KPAR:-1}; [ "$KPAR" -lt 1 ] && KPAR=1
    NP=$SLURM_NTASKS
    [ -n "$NIONS" ] && [ "$NIONS" -lt "$NP" ] && NP=$NIONS
    NP=$(( NP / KPAR * KPAR ))
    [ "$NP" -lt "$KPAR" ] && NP=$KPAR
    echo "== $JID step $STEP_NUM: running VASP with $NP ranks ($NIONS ions, KPAR=$KPAR) in $STEP_DIR ($(date))"
    build_potcar "$STEP_DIR"
    ( cd "$STEP_DIR" && mpirun -np "$NP" "$VASP_BIN" )
    mkdir -p "$STEP_DIR/results"
    cp "$STEP_DIR"/OUTCAR "$STEP_DIR"/CONTCAR "$STEP_DIR/results/"
else
    echo "== $JID step $STEP_NUM: results/OUTCAR already present, skipping VASP"
fi

echo "== $JID: recording step $STEP_NUM"
python dft_prep.py record "$JID" "$STEP_NUM"

echo "== $JID: generating next step"
if python dft_prep.py next "$JID"; then
    NEW_STEP=$(ls -d "$SUP_DIR"/step_* | sort -t_ -k2 -n | tail -1)
    if [ "$NEW_STEP" = "$STEP_DIR" ]; then
        echo "== $JID: next reported success but produced no new step, stopping to avoid a loop"
        python dft_prep.py voltage "$JID" || true
        exit 0
    fi
    build_potcar "$NEW_STEP"
    echo "== $JID: resubmitting chain for $(basename "$NEW_STEP")"
    sbatch --job-name="bm_$JID" "$0" "$JID"
else
    echo "== $JID: no further step (delithiation complete), computing voltage curve"
    python dft_prep.py voltage "$JID" || true
fi
echo "== $JID: chain segment done ($(date))"
