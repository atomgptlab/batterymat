"""Task 2 statistics: tier-1 vs tier-2 agreement over the Li_ screening pool.

Reads tier1_pool_predictions.csv, writes pool_stats.json and pool_parity.png.
v_tier2 = Li_min.csv avg_voltage (screening-pipeline ALIGNN-FF labels);
v_label = summary.csv avg_voltage_V (the ALIGNN-FF label set the tier-1 model
was trained on). These are two different ALIGNN-FF runs, so all three
pairwise agreements are reported.
"""
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 15,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#52514e"
LO, HI = 3.0, 4.5


def agree(a, b):
    m = a.notna() & b.notna()
    if m.sum() < 3:
        return {"n": int(m.sum())}
    a, b = a[m].values, b[m].values
    d = a - b
    ia, ib = (a >= LO) & (a <= HI), (b >= LO) & (b <= HI)
    return {"n": int(m.sum()), "spearman": float(spearmanr(a, b).correlation), "pearson": float(pearsonr(a, b)[0]),
            "mae_V": float(np.mean(np.abs(d))), "rmse_V": float(np.sqrt(np.mean(d ** 2))),
            "mean_diff_V(first-second)": float(np.mean(d)), "median_diff_V": float(np.median(d)),
            "window_agreement": float(np.mean(ia == ib)), "n_first_in_window": int(ia.sum()), "n_second_in_window": int(ib.sum()),
            "n_both_in_window": int((ia & ib).sum())}


def main():
    df = pd.read_csv(os.path.join(HERE, "tier1_pool_predictions.csv"))
    out = {"pool_size": int(len(df)), "n_pred": int(df["v_tier1"].notna().sum()),
           "partition_counts": df["partition"].value_counts().to_dict(),
           "n_with_v_label": int(df["v_label"].notna().sum())}
    notrain = df[df["partition"] != "train"]
    test = df[df["partition"] == "test"]
    none = df[df["partition"] == "none"]
    out["tier1_vs_tier2_all"] = agree(df["v_tier1"], df["v_tier2"])
    out["tier1_vs_tier2_not_in_train"] = agree(notrain["v_tier1"], notrain["v_tier2"])
    out["tier1_vs_tier2_test_only"] = agree(test["v_tier1"], test["v_tier2"])
    out["tier1_vs_tier2_not_in_training_set_at_all(none)"] = agree(none["v_tier1"], none["v_tier2"])
    out["tier1_vs_label_all"] = agree(df["v_tier1"], df["v_label"])
    out["tier1_vs_label_not_in_train"] = agree(notrain["v_tier1"], notrain["v_label"])
    out["tier1_vs_label_test_only"] = agree(test["v_tier1"], test["v_label"])
    out["label_vs_tier2_all"] = agree(df["v_label"], df["v_tier2"])
    out["label_vs_tier2_not_in_train"] = agree(notrain["v_label"], notrain["v_tier2"])
    # in-window-by-tier-2 subset (screening relevant band)
    band = df[(df["v_tier2"] >= LO) & (df["v_tier2"] <= HI)]
    out["tier1_vs_tier2_within_tier2_window"] = agree(band["v_tier1"], band["v_tier2"])
    json.dump(out, open(os.path.join(HERE, "pool_stats.json"), "w"), indent=2)
    print(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    for ax, (xcol, xl, key) in zip(axes, [("v_tier2", "Tier-2 ALIGNN-FF, Li_min.csv (V)", "tier1_vs_tier2_all"),
                                          ("v_label", "ALIGNN-FF training label, summary.csv (V)", "tier1_vs_label_all")]):
        m = df[xcol].notna() & df["v_tier1"].notna()
        tr = m & (df["partition"] == "train")
        ax.scatter(df.loc[tr, xcol], df.loc[tr, "v_tier1"], s=6, color=BLUE, alpha=0.35, edgecolors="none", label=f"train ({tr.sum()})")
        nt = m & (df["partition"] != "train")
        ax.scatter(df.loc[nt, xcol], df.loc[nt, "v_tier1"], s=8, color=ORANGE, alpha=0.6, edgecolors="none", label=f"val/test/none ({nt.sum()})")
        lim = [-3, 8]
        ax.plot(lim, lim, color=GRAY, lw=1.2); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.axvspan(LO, HI, color=ORANGE, alpha=0.08); ax.axhspan(LO, HI, color=ORANGE, alpha=0.08)
        s = out[key]
        ax.set_title(f"n={s['n']}: ρ={s['spearman']:.3f}, r={s['pearson']:.3f}, MAE={s['mae_V']:.3f} V")
        ax.set_xlabel(xl); ax.set_ylabel("Tier-1 prediction (V)")
        ax.legend(frameon=False, loc="upper left", markerscale=3)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "pool_parity.png"), dpi=200)


if __name__ == "__main__":
    main()
