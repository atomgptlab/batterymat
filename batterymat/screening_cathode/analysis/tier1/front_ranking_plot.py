"""Task 3 figure: tier-1 vs tier-2 rank for common survivors + voltages.
Writes front_ranking.png (reads tier1_front_ranking.csv, front_ranking_comparison.json)."""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 15,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13})
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#52514e"
t1 = pd.read_csv(os.path.join(HERE, "tier1_front_ranking.csv"))
cmp = json.load(open(os.path.join(HERE, "front_ranking_comparison.json")))
ranked = pd.read_csv(os.path.abspath(os.path.join(HERE, "..", "..", "cathode_candidates_ranked.csv")))
pool = pd.read_csv(os.path.join(HERE, "tier1_pool_predictions.csv")).set_index("jid")

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
ax = axes[0]
c = t1.dropna(subset=["rank_tier2"])
ax.scatter(c["rank_tier2"], c["rank"], s=16, color=BLUE, alpha=0.7, edgecolors="none", label=f"common survivors ({len(c)})")
lim = [0, max(len(t1), len(ranked)) + 10]
ax.plot(lim, lim, color=GRAY, lw=1.2, ls="--")
for j, lab in [("JVASP-2017", "LCO"), ("JVASP-117295", "Li3FeO3"), ("JVASP-96563", "Li2MgMn3O8"), ("JVASP-116849", "LiV2F7")]:
    r = c[c["jid"] == j]
    if len(r):
        ax.scatter(r["rank_tier2"], r["rank"], s=90, color=ORANGE, edgecolors="black", zorder=3)
        ax.annotate(lab, (r["rank_tier2"].values[0], r["rank"].values[0]), xytext=(8, -4), textcoords="offset points", fontsize=13)
ax.set_xlabel("Rank, tier-2 ALIGNN-FF front (n = %d)" % len(ranked))
ax.set_ylabel("Rank, tier-1 ALIGNN front (n = %d)" % len(t1))
ax.set_title(f"Spearman ρ = {cmp['spearman_rank_common']:.3f}; top-12 overlap {cmp['overlap_top12']['n_common']}/12, top-50 {cmp['overlap_top50']['n_common']}/50")
ax.legend(frameon=False, loc="lower right")
ax = axes[1]
surv = set(t1["jid"]) | set(ranked["jid"])
s = pool.loc[sorted(surv)]
both = s.index.isin(t1["jid"]) & s.index.isin(ranked["jid"])
only1 = s.index.isin(t1["jid"]) & ~s.index.isin(ranked["jid"])
only2 = ~s.index.isin(t1["jid"]) & s.index.isin(ranked["jid"])
ax.scatter(s.loc[both, "v_tier2"], s.loc[both, "v_tier1"], s=16, color=BLUE, alpha=0.6, edgecolors="none", label=f"both fronts ({both.sum()})")
ax.scatter(s.loc[only2, "v_tier2"], s.loc[only2, "v_tier1"], s=22, color=ORANGE, alpha=0.8, edgecolors="none", label=f"tier-2 only ({only2.sum()})")
ax.scatter(s.loc[only1, "v_tier2"], s.loc[only1, "v_tier1"], s=22, color="#1baf7a", alpha=0.8, edgecolors="none", label=f"tier-1 only ({only1.sum()})")
ax.axhspan(3.0, 4.5, color=GRAY, alpha=0.08); ax.axvspan(3.0, 4.5, color=GRAY, alpha=0.08)
ax.plot([1.5, 5.5], [1.5, 5.5], color=GRAY, lw=1.2, ls="--")
ax.set_xlim(1.5, 5.5); ax.set_ylim(1.5, 5.5)
ax.set_xlabel("Tier-2 ALIGNN-FF average voltage (V)"); ax.set_ylabel("Tier-1 ALIGNN average voltage (V)")
ax.set_title("Survivors of either front; shaded = 3.0–4.5 V window")
ax.legend(frameon=False, loc="upper left")
for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "front_ranking.png"), dpi=200)
print("saved")
