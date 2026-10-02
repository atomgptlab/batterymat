"""Task 5 (R1.m5): training vs validation loss curves from history_*.json.

Each history entry is [total_loss, graphwise_loss, 0, 0, 0, 0] per epoch
(MSE criterion, summed over batches -> train and val totals are on different
scales because the train split has 8x the batches; both are also plotted
normalised per batch). Output: loss_curves.png, loss_curves.json
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = "/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/voltage"
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 16,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})
BLUE, ORANGE = "#2a78d6", "#eb6834"

tr = json.load(open(os.path.join(MODEL_DIR, "history_train.json")))
va = json.load(open(os.path.join(MODEL_DIR, "history_val.json")))
ids = json.load(open(os.path.join(MODEL_DIR, "ids_train_val_test.json")))
bs = 32
# number of batches per epoch = number of batched entries in the archived
# Train_results.json / Val_results.json (drop_last=True: 190 and 23)
tr_res = json.load(open(os.path.join(MODEL_DIR, "Train_results.json")))
va_res = json.load(open(os.path.join(MODEL_DIR, "Val_results.json")))
n_tr_batches, n_va_batches = len(tr_res), len(va_res)
def _mae(res):
    t = np.concatenate([np.array(r["target_out"]).ravel() for r in res])
    p = np.concatenate([np.array(r["pred_out"]).ravel() for r in res])
    return float(np.mean(np.abs(p - t))), float(np.sqrt(np.mean((p - t) ** 2))), int(len(t))
tr_mae, tr_rmse, n_tr = _mae(tr_res); va_mae, va_rmse, n_va = _mae(va_res)
trl = np.array([e[0] for e in tr]); val = np.array([e[0] for e in va])
ep = np.arange(1, len(trl) + 1)
trn = trl / n_tr_batches; van = val / n_va_batches  # mean MSE per batch
best = int(np.argmin(val)) + 1

fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axes[0]
ax.plot(ep, trl, color=BLUE, lw=2, label="train (summed over batches)")
ax.plot(ep, val, color=ORANGE, lw=2, label="validation (summed over batches)")
ax.set_yscale("log"); ax.set_xlabel("Epoch"); ax.set_ylabel("MSE loss, summed over batches")
ax.set_title("Summed over batches (as logged)")
ax.legend(frameon=False)
ax = axes[1]
ax.plot(ep, trn, color=BLUE, lw=2, label=f"train ({n_tr_batches} batches/epoch)")
ax.plot(ep, van, color=ORANGE, lw=2, label=f"validation ({n_va_batches} batches/epoch)")
ax.axvline(best, color="#52514e", lw=1.2, ls="--")
ax.set_yscale("log"); ax.set_xlabel("Epoch"); ax.set_ylabel("Mean MSE per batch (V²)")
ax.set_title(f"Per-batch; best val at epoch {best}")
ax.legend(frameon=False)
for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "loss_curves.png"), dpi=200)

summary = {"epochs": len(trl), "batch_size": bs, "n_train_batches": n_tr_batches, "n_val_batches": n_va_batches,
           "best_val_epoch": best, "train_loss_summed": {str(e): float(trl[e-1]) for e in (1, 50, 100, 150, 200, 250)},
           "val_loss_summed": {str(e): float(val[e-1]) for e in (1, 50, 100, 150, 200, 250)},
           "train_mse_per_batch_final": float(trn[-1]), "val_mse_per_batch_final": float(van[-1]),
           "train_rmse_per_batch_final_V": float(np.sqrt(trn[-1])), "val_rmse_per_batch_final_V": float(np.sqrt(van[-1])),
           "val_over_train_per_batch_final": float(van[-1] / trn[-1]),
           "train_results_json": {"n": n_tr, "mae_V": tr_mae, "rmse_V": tr_rmse},
           "val_results_json": {"n": n_va, "mae_V": va_mae, "rmse_V": va_rmse}}
json.dump(summary, open(os.path.join(HERE, "loss_curves.json"), "w"), indent=2)
print(json.dumps(summary, indent=1))
