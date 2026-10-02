"""Tasks 1 + 4: held-out test metrics, full reproduction check, screening
confusion matrix (3.0 <= V <= 4.5 window), residual distribution + bootstrap.

Outputs: test_metrics.json, residuals.png, test_reproduction.csv
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict import ARCHIVE, MODEL_DIR, Tier1Predictor  # noqa: E402

plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 16,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#52514e"
LO, HI = 3.0, 4.5


def main():
    test = json.load(open(os.path.join(MODEL_DIR, "Test_results.json")))
    y = np.array([r["target_out"][0] for r in test])
    yhat = np.array([r["pred_out"] for r in test])
    ids = [r["id"] for r in test]
    n = len(y)

    # --- full reproduction of the archived predictions (all 761) ---
    idprop = {d["jid"]: d for d in json.load(open(os.path.join(ARCHIVE, "id_prop.json")))}
    p = Tier1Predictor()
    rec = np.array([p.predict_one(idprop[j]["atoms"]) for j in ids])
    repro = pd.DataFrame({"jid": ids, "label": y, "pred_archived": yhat, "pred_recomputed": rec,
                          "abs_diff": np.abs(rec - yhat)})
    repro.to_csv(os.path.join(HERE, "test_reproduction.csv"), index=False)

    # --- metrics (task 1) ---
    res = yhat - y
    mae = float(np.mean(np.abs(res)))
    rmse = float(np.sqrt(np.mean(res ** 2)))
    r2 = float(1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2))
    mad = float(np.mean(np.abs(y - y.mean())))
    # MAD quoted in the archived `mad` file is computed on the whole dataset
    ip_targets = np.array([d["target"] for d in idprop.values()])
    mad_all = float(np.mean(np.abs(ip_targets - ip_targets.mean())))
    pear = float(np.corrcoef(y, yhat)[0, 1])
    from scipy.stats import spearmanr
    spear = float(spearmanr(y, yhat).correlation)

    # --- confusion matrix on the 3.0-4.5 V window (task 4 / R2.M8) ---
    in_lab = (y >= LO) & (y <= HI)
    in_pred = (yhat >= LO) & (yhat <= HI)
    tp = int(np.sum(in_lab & in_pred)); fp = int(np.sum(~in_lab & in_pred))
    fn = int(np.sum(in_lab & ~in_pred)); tn = int(np.sum(~in_lab & ~in_pred))
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec_ = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec_ / (prec + rec_) if prec + rec_ else float("nan")
    acc = (tp + tn) / n
    # false negatives / positives near the window edges
    fn_ids = [ids[i] for i in np.where(in_lab & ~in_pred)[0]]
    fp_ids = [ids[i] for i in np.where(~in_lab & in_pred)[0]]

    # --- residual distribution + bootstrap ---
    rng = np.random.default_rng(0)
    boots = np.array([np.mean(np.abs(res[rng.integers(0, n, n)])) for _ in range(10000)])
    ci = np.percentile(boots, [2.5, 97.5])
    absres = np.abs(res)
    out = {
        "n_test": n,
        "reproduction": {"n": n, "max_abs_diff_V": float(repro["abs_diff"].max()),
                         "mean_abs_diff_V": float(repro["abs_diff"].mean()),
                         "n_above_1e-3": int((repro["abs_diff"] > 1e-3).sum())},
        "mae_V": mae, "rmse_V": rmse, "r2": r2, "pearson_r": pear, "spearman_rho": spear,
        "mad_test_V": mad, "mad_all_data_V": mad_all, "mad_archived_file_V": 1.2443734639220472,
        "mae_over_mad_test": mae / mad,
        "residual_mean_V": float(res.mean()), "residual_std_V": float(res.std(ddof=1)),
        "residual_median_V": float(np.median(res)),
        "abs_residual_p50_V": float(np.percentile(absres, 50)),
        "abs_residual_p90_V": float(np.percentile(absres, 90)),
        "abs_residual_p95_V": float(np.percentile(absres, 95)),
        "abs_residual_p99_V": float(np.percentile(absres, 99)),
        "abs_residual_max_V": float(absres.max()),
        "signed_residual_2.5_97.5_pct_V": [float(v) for v in np.percentile(res, [2.5, 97.5])],
        "signed_residual_0.5_99.5_pct_V": [float(v) for v in np.percentile(res, [0.5, 99.5])],
        "frac_abs_residual_le_0.25V": float(np.mean(absres <= 0.25)),
        "frac_abs_residual_le_0.5V": float(np.mean(absres <= 0.5)),
        "bootstrap_mae_95ci_V": [float(ci[0]), float(ci[1])], "bootstrap_resamples": 10000, "bootstrap_seed": 0,
        "window_V": [LO, HI],
        "confusion": {"TP": tp, "FP": fp, "FN": fn, "TN": tn, "precision": prec, "recall": rec_,
                      "f1": f1, "accuracy": acc, "n_in_window_label": int(in_lab.sum()),
                      "n_in_window_pred": int(in_pred.sum()), "FN_ids": fn_ids, "FP_ids": fp_ids},
    }
    json.dump(out, open(os.path.join(HERE, "test_metrics.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "confusion"}, indent=1))
    print("confusion", {k: v for k, v in out["confusion"].items() if not k.endswith("ids")})

    # --- figure: residual histogram + parity ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    ax = axes[0]
    ax.hist(res, bins=60, color=BLUE, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color=GRAY, lw=1.5)
    for q in np.percentile(res, [2.5, 97.5]):
        ax.axvline(q, color=ORANGE, lw=2, ls="--")
    ax.set_xlabel("Residual  V$_{pred}$ $-$ V$_{label}$ (V)")
    ax.set_ylabel("Count")
    ax.set_title(f"Held-out test (n={n}): MAE {mae:.3f} V, σ {res.std(ddof=1):.3f} V")
    ax.text(0.02, 0.97, f"dashed: 2.5/97.5 pct\n[{np.percentile(res,2.5):+.2f}, {np.percentile(res,97.5):+.2f}] V",
            transform=ax.transAxes, va="top", fontsize=13, color=ORANGE)
    ax = axes[1]
    ax.scatter(y, yhat, s=14, color=BLUE, alpha=0.6, edgecolors="none")
    lim = [min(y.min(), yhat.min()) - 0.3, max(y.max(), yhat.max()) + 0.3]
    ax.plot(lim, lim, color=GRAY, lw=1.5)
    ax.axvspan(LO, HI, color=ORANGE, alpha=0.10)
    ax.axhspan(LO, HI, color=ORANGE, alpha=0.10)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("ALIGNN-FF label (V)"); ax.set_ylabel("Tier-1 prediction (V)")
    ax.set_title(f"R² = {r2:.3f};  3.0–4.5 V window: P {prec:.2f}, R {rec_:.2f}")
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "residuals.png"), dpi=200)
    print("saved residuals.png")


if __name__ == "__main__":
    main()
