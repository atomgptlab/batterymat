"""Item 9: ALIGNN tier-1 parity plot from the archived Test_results.json.

Predicted vs target (ALIGNN-FF) average voltage on the 761 held-out test
structures; dashed y=x, dotted dataset-mean baseline (the predictor whose MAE
is the MAD), MAE and R2 in the panel. Writes alignn_parity.png / .pdf.
Numbers match test_metrics.json (same file, same formulas).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = "/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/voltage"
plt.rcParams.update({"font.size": 15, "axes.labelsize": 17, "axes.titlesize": 16,
                     "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 14,
                     "pdf.fonttype": 42})
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#52514e"

test = json.load(open(os.path.join(MODEL_DIR, "Test_results.json")))
y = np.array([r["target_out"][0] for r in test]); yhat = np.array([r["pred_out"] for r in test])
res = yhat - y
mae = np.mean(np.abs(res)); r2 = 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)
mad = np.mean(np.abs(y - y.mean()))
lim = [np.floor(min(y.min(), yhat.min())) - 0.5, np.ceil(max(y.max(), yhat.max())) + 0.5]

fig, ax = plt.subplots(figsize=(6.2, 6.2))
ax.plot(lim, lim, ls="--", color=GRAY, lw=1.6, label="y = x", zorder=1)
ax.axhline(y.mean(), ls=":", color=ORANGE, lw=2.0, label=f"dataset mean ({y.mean():.2f} V)", zorder=1)
ax.scatter(y, yhat, s=18, color=BLUE, alpha=0.6, edgecolors="none", label=f"held-out test (n = {len(y)})", zorder=2)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
ax.set_xlabel("ALIGNN-FF average voltage (V)")
ax.set_ylabel("ALIGNN predicted average voltage (V)")
ax.text(0.04, 0.96, f"MAE = {mae:.3f} V\nR$^2$ = {r2:.3f}\nMAD (test) = {mad:.3f} V", transform=ax.transAxes,
        va="top", ha="left", fontsize=15, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=GRAY, lw=0.8))
ax.legend(loc="lower right", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "alignn_parity.png"), dpi=300)
fig.savefig(os.path.join(HERE, "alignn_parity.pdf"))
print(f"n={len(y)} MAE={mae:.4f} R2={r2:.4f} MAD={mad:.4f} mean={y.mean():.4f}")
