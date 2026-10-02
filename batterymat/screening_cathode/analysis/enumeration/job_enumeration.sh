#!/bin/bash
#SBATCH --job-name=bm_enum
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G          # without this atomgptlab allocates the whole node (DefMemPerNode=UNLIMITED)
#SBATCH --time=24:00:00
#SBATCH --output=enum_%j.out
#
# Exhaustive Li-vacancy enumeration with ALIGNN-FF (reviewer R1.M6), level (b):
# relax all 2^8 configurations of LiCoO2 (JVASP-2017, 2x2x2) and Li2MgMn3O8
# (JVASP-96563, 1x1x1) and rebuild the summary. Level (a) single points are
# recomputed here for clean provenance; the Mac results_*.json are kept as
# results_*_mac.json (their energies match atomgptlab to <0.05 meV, but the Mac
# alignn 2025.4.1 forces are wrong, so relaxations must run here).
#
# Submit from this directory on the AtomGPT-lab cluster:
#   sbatch job_enumeration.sh
# Requirements on the cluster:
#   - conda env `alignn` with alignn + jarvis-tools + ase (the env that ran the
#     August prospective campaign); ALIGNN-FF alignnff_wt01 weights available
#   - JARVIS dft_3d pickle containing JVASP-913, JVASP-2017, JVASP-96563; set CACHE
#   - POSCAR-JVASP-96563.vasp (reference cell) next to this script or set REF96563
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
CONDA_SH="${CONDA_SH:-/data/$USER/miniforge3/etc/profile.d/conda.sh}"   # atomgptlab
[ -f "$CONDA_SH" ] || CONDA_SH="$HOME/miniforge3/etc/profile.d/conda.sh"   # skipjack
source "$CONDA_SH"
conda activate alignn
# Skipjack runs ALIGNN 2.0 from a source checkout; atomgptlab has alignn installed in the env
[ -d "${ALIGNN_REPO:-$HOME/alignn/repo}/alignn" ] && export PYTHONPATH="${ALIGNN_REPO:-$HOME/alignn/repo}${PYTHONPATH:+:$PYTHONPATH}"
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
export CACHE="${CACHE:-/data/$USER/Li_250/jarvis_dft3d_cache.pkl}"
export REF96563="${REF96563:-./POSCAR-JVASP-96563.vasp}"
W=${SLURM_CPUS_PER_TASK:-32}

# Sanity check of the force field before spending CPU time: forces on a
# DFT-relaxed cell must be small. Abort if fmax > 2 eV/A.
python - <<'PY'
import os, pickle, numpy as np
from jarvis.core.atoms import Atoms
from alignn.ff.ff import AlignnAtomwiseCalculator, wt01_path
d = pickle.load(open(os.environ["CACHE"], "rb"))
x = [e for e in d if e["jid"] == "JVASP-2017"][0]
a = Atoms.from_dict(x["atoms"]).ase_converter()
a.calc = AlignnAtomwiseCalculator(path=wt01_path(), stress_wt=0.3)
f = np.abs(a.get_forces()).max()
print(f"ALIGNN-FF fmax on DFT-relaxed LiCoO2: {f:.3f} eV/A")
raise SystemExit(0 if f < 2.0 else 1)
PY

for j in JVASP-2017 JVASP-96563; do [ -f results_$j.json ] && [ ! -f results_${j}_mac.json ] && mv results_$j.json results_${j}_mac.json || true; done
python -u enumerate.py JVASP-2017 --resume --relax all --workers "$W" --cache "$CACHE"
python -u enumerate.py JVASP-96563 --resume --relax all --workers "$W" --cache "$CACHE" --ref-poscar "$REF96563"
python enumerate.py --summary
