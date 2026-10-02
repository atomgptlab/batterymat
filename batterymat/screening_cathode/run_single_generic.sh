#!/bin/bash
#SBATCH --job-name=bmsingle
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=single_%x_%j.out
#SBATCH --error=single_%x_%j.err

# =============================================================================
# One VASP relaxation in one directory (no chain). Used for
#   campaign C  : dft_inputs/spotcheck/<JID>/{lithiated,delithiated}
#   references  : dft_inputs/ref-Na/Na_pv_{PBE,optB88vdW}, dft_inputs/ref-Mg/Mg_pv_{PBE,optB88vdW}
# Builds POTCAR from POTCAR_spec, runs vasp_std, copies OUTCAR/CONTCAR/OSZICAR
# to results/. Edit the block below and the #SBATCH header for your cluster.
#
# Usage (from $WORKDIR):
#   sbatch --job-name=sc_JVASP-156413_lith run_single_generic.sh dft_inputs/spotcheck/JVASP-156413/lithiated
#   for d in dft_inputs/spotcheck/*/*/; do n=$(echo $d | awk -F/ '{print $3"_"$4}');
#       sbatch --job-name="sc_$n" run_single_generic.sh "$d"; done
#   for d in dft_inputs/ref-*/*/; do sbatch --job-name="ref_$(basename $d)" run_single_generic.sh "$d"; done
# =============================================================================

# ----------------------------- EDIT ME ---------------------------------------
WORKDIR="${BMAT_WORKDIR:-/data/$USER/revision}"
VASP_MODULE="${BMAT_VASP_MODULE:-}"
VASP_BIN="${BMAT_VASP_BIN:-}"
POTCAR_DIR="${BMAT_POTCAR_DIR:-${VASP_PP_PATH:-}}"
VDW_KERNEL="${BMAT_VDW_KERNEL:-}"
# -----------------------------------------------------------------------------

set -e
set -o pipefail
RUN_DIR="$1"
[ -z "$RUN_DIR" ] && { echo "usage: sbatch run_single_generic.sh <run directory>"; exit 1; }
cd "$WORKDIR"
[ -d "$RUN_DIR" ] || { echo "FATAL: $RUN_DIR not found under $WORKDIR"; exit 1; }

export OMP_NUM_THREADS=1
ulimit -s unlimited
[ -n "$VASP_MODULE" ] && module load "$VASP_MODULE"
[ -z "$VASP_BIN" ] && command -v vasp_std >/dev/null 2>&1 && VASP_BIN=$(command -v vasp_std)
[ -z "$VASP_BIN" ] && { echo "FATAL: vasp_std not found. Set VASP_BIN."; exit 1; }
[ -z "$POTCAR_DIR" ] || [ ! -d "$POTCAR_DIR" ] && { echo "FATAL: POTCAR_DIR not set/found."; exit 1; }

if [ -f "$RUN_DIR/results/OUTCAR" ]; then
    echo "== $RUN_DIR: results/OUTCAR already present, nothing to do"; exit 0
fi

if [ ! -s "$RUN_DIR/POTCAR" ]; then
    : > "$RUN_DIR/POTCAR"
    while read -r el; do
        [ -z "$el" ] && continue
        cat "$POTCAR_DIR/$el/POTCAR" >> "$RUN_DIR/POTCAR"
    done < "$RUN_DIR/POTCAR_spec"
fi
if grep -qi "LUSE_VDW" "$RUN_DIR/INCAR" && [ -n "$VDW_KERNEL" ] && [ ! -f "$RUN_DIR/vdw_kernel.bindat" ]; then
    cp "$VDW_KERNEL" "$RUN_DIR/vdw_kernel.bindat"
fi

# Rank count: at most one rank per ion, multiple of KPAR (small cells otherwise abort in M_divide)
NIONS=$(awk 'NR==7{for(i=1;i<=NF;i++)s+=$i; print s}' "$RUN_DIR/POSCAR")
KPAR=$(awk -F= 'toupper($0) ~ /^ *KPAR/ {print $2+0}' "$RUN_DIR/INCAR" | head -1)
KPAR=${KPAR:-1}; [ "$KPAR" -lt 1 ] && KPAR=1
NP=$SLURM_NTASKS
[ -n "$NIONS" ] && [ "$NIONS" -lt "$NP" ] && NP=$NIONS
NP=$(( NP / KPAR * KPAR )); [ "$NP" -lt "$KPAR" ] && NP=$KPAR
# NCORE must divide the ranks per k-point group, or VASP aborts in M_divide on small cells.
NC=$(awk -F= 'toupper($0) ~ /^ *NCORE/ {print $2+0}' "$RUN_DIR/INCAR" | head -1)
if [ -n "$NC" ] && [ "$NC" -gt 1 ]; then
    PER=$(( NP / KPAR ))
    while [ $(( PER % NC )) -ne 0 ]; do NC=$(( NC - 1 )); done
    sed -i "s/^ *NCORE *=.*/NCORE = $NC/" "$RUN_DIR/INCAR"
fi
echo "== $RUN_DIR: VASP with $NP ranks ($NIONS ions, KPAR=$KPAR) ($(date))"
( cd "$RUN_DIR" && mpirun -np "$NP" "$VASP_BIN" )
mkdir -p "$RUN_DIR/results"
cp "$RUN_DIR"/OUTCAR "$RUN_DIR"/CONTCAR "$RUN_DIR"/OSZICAR "$RUN_DIR/results/"
rm -f "$RUN_DIR"/WAVECAR "$RUN_DIR"/CHG
grep -q "reached required accuracy" "$RUN_DIR/results/OUTCAR" || echo "NOT_CONVERGED: $RUN_DIR hit NSW without reaching EDIFFG; energy recorded anyway, inspect before use"
grep "free  energy   TOTEN" "$RUN_DIR/results/OUTCAR" | tail -1
echo "== $RUN_DIR: done ($(date))"
