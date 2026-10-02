#!/bin/bash
#SBATCH --job-name=bm_oodcv
#SBATCH --partition=gpu
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=oodcv_%j.out
#
# Grouped-split generalisation study for the tier-1 ALIGNN regressor (R2.M7).
# 1. make_grouped_splits.py writes fold directories (leave-one-redox-metal-out,
#    leave-one-anion-class-out, and 5 random folds) from the archived id_prop.json.
# 2. Each fold is trained with the archived config (Li_250/config.json) but with
#    EPOCHS epochs (default 100 to fit the budget; state this in the paper) using
#    ALIGNN's train_alignn with explicit train/val/test files.
# 3. collect_grouped_results.py gathers per-fold test MAE / R2 into oodcv_results.csv.
#
# Set ARCHIVE to the copied Li_250 directory on the cluster (needs id_prop.json and
# config.json). Copy from the Mac:
#   rsync -av ~/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/{id_prop.json,config.json} <cluster>:/data/$USER/Li_250/
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ "$(uname -m)" = aarch64 ]; then
    # atomgptlab GPU nodes (NVIDIA GB10, aarch64): same env/bundle as the ffF30 training jobs
    PY=/data/$USER/miniforge3-aarch64/envs/alignn/bin/python
    export PYTHONPATH="${ALIGNN_BUNDLE:-/data/$USER/bundles/alignn_12ae44e}"
else
    CONDA_SH="${CONDA_SH:-/data/$USER/miniforge3/etc/profile.d/conda.sh}"   # atomgptlab
    [ -f "$CONDA_SH" ] || CONDA_SH="$HOME/miniforge3/etc/profile.d/conda.sh"   # skipjack
    source "$CONDA_SH"
    conda activate alignn
    # Skipjack runs ALIGNN 2.0 from a source checkout; atomgptlab has alignn installed in the env
    [ -d "${ALIGNN_REPO:-$HOME/alignn/repo}/alignn" ] && export PYTHONPATH="${ALIGNN_REPO:-$HOME/alignn/repo}${PYTHONPATH:+:$PYTHONPATH}"
    PY=python
fi
ARCHIVE="${ARCHIVE:-/data/$USER/Li_250}"
EPOCHS="${EPOCHS:-100}"

"$PY" make_grouped_splits.py --id-prop "$ARCHIVE/id_prop.json" --out folds
for d in folds/*/; do
  name=$(basename "$d")
  [ -f "$d/out/prediction_results_test_set.csv" ] && { echo "skip $name (done)"; continue; }
  "$PY" - "$ARCHIVE/config.json" "$d" "$EPOCHS" <<'PY'
import json, sys
d = sys.argv[2].rstrip("/")
cfg = json.load(open(sys.argv[1])); info = json.load(open(d + "/fold_info.json"))
cfg["epochs"] = int(sys.argv[3]); cfg["output_dir"] = d + "/out"
cfg["n_train"], cfg["n_val"], cfg["n_test"] = info["n_train"], info["n_val"], info["n_test"]
cfg["keep_data_order"] = True
if not cfg["model"].get("calculate_gradient", False):
    # gradient loss is inactive without calculate_gradient; zero its weight so
    # train.py writes prediction_results_test_set.csv (it requires gradwise_weight == 0)
    cfg["model"]["gradwise_weight"] = 0.0
import importlib.util
if importlib.util.find_spec("dgl") is None:
    # DGL-free ALIGNN (atomgptlab GB10 bundle): pure-torch graph and model with the
    # same layer sizes; cutoff 8 A / 12 neighbours as in Li_250, and the three-body
    # cutoff widened from the 3.5 A default to the full cutoff so the line graph
    # spans all neighbour pairs as in the original k-nearest DGL graph.
    cfg["neighbor_strategy"] = "pure_torch"
    cfg["model"]["name"] = "alignn_atomwise_pure"
    cfg["three_body_cutoff"] = cfg.get("cutoff", 8.0)
json.dump(cfg, open(d + "/config.json", "w"), indent=1)
PY
  D="$(cd "$d" && pwd)"; mkdir -p "$D/work"
  (cd "$D/work" && "$PY" -m alignn.train_alignn --root_dir "$D" --config_name "$D/config.json" \
      --target_key target --id_key jid --output_dir "$D/out") \
      || echo "fold $name failed"
done
"$PY" collect_grouped_results.py folds > oodcv_results.csv
