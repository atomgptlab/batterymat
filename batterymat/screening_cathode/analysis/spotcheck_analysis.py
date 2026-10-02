"""Statistics for the campaign-C DFT spot-check of the tier-1 ALIGNN labels.

Reads spotcheck_results.csv (written by spotcheck_collect.py), marks a pair as
unconverged when either endpoint OUTCAR lacks "reached required accuracy", joins
the tier-1 regressor prediction for the same held-out structures
(tier1/test_reproduction.csv), and reports, for the DFT voltage versus the
ALIGNN-FF label and versus the regressor prediction:

  - MAE, mean error (bias), RMSE, and MAE after removing the mean bias
  - Spearman rank correlation
  - agreement of the screen's in-window call (3.0 <= V <= 4.5 V)

for (i) all converged pairs and (ii) the converged pairs whose formula contains
a redox-active transition metal (the screen's `REDOX_METALS` criterion, which
all 682 ranked candidates satisfy). Writes spotcheck_stats.json.

Usage (from screening_cathode/analysis/):  python spotcheck_analysis.py
"""
import csv
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPOT = HERE.parent / "dft_inputs" / "spotcheck"
RESTART = HERE.parent / "dft_inputs" / "extra_2026-10-01" / "spotcheck-restart"
REDOX_METALS = {"Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                "Nb", "Mo", "Ru", "Rh", "W"}   # same set as screen_cathode.py
V_WINDOW = (3.0, 4.5)


def converged(run_dir: Path) -> bool:
    for p in (run_dir / "results" / "OUTCAR", run_dir / "OUTCAR"):
        if p.exists():
            return "reached required accuracy" in p.read_text(errors="replace")
    return False


def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def stats(dft, ref):
    dft, ref = np.asarray(dft), np.asarray(ref)
    e = dft - ref
    inw = lambda v: (v >= V_WINDOW[0]) & (v <= V_WINDOW[1])
    return {
        "n": int(len(e)),
        "MAE_V": round(float(np.mean(np.abs(e))), 3),
        "ME_V": round(float(np.mean(e)), 3),
        "RMSE_V": round(float(np.sqrt(np.mean(e ** 2))), 3),
        "MAE_bias_removed_V": round(float(np.mean(np.abs(e - e.mean()))), 3),
        "spearman": round(spearman(dft, ref), 3),
        "window_agreement": f"{int(np.sum(inw(dft) == inw(ref)))}/{len(e)}",
    }


def main():
    rows = list(csv.DictReader(open(HERE / "spotcheck_results.csv")))
    pred = {r["jid"]: float(r["pred_recomputed"])
            for r in csv.DictReader(open(HERE / "tier1" / "test_reproduction.csv"))}
    table = []
    for r in rows:
        if not r["dft_V"]:
            continue
        def side(name):  # restart directory (2026-10-01) replaces an endpoint that hit NSW
            rd = RESTART / r["jid"] / name
            return rd if (rd / "results" / "OUTCAR").exists() else SPOT / r["jid"] / name
        conv = converged(side("lithiated")) and converged(side("delithiated"))
        els = set(re.findall(r"[A-Z][a-z]?", r["formula"]))
        table.append(dict(jid=r["jid"], formula=r["formula"], functional=r["functional"],
                          label=float(r["label_V"]), pred=pred[r["jid"]], dft=float(r["dft_V"]),
                          converged=conv, redox=bool(els & REDOX_METALS)))

    out = {"excluded_unconverged": [f"{t['jid']} {t['formula']}" for t in table if not t["converged"]]}
    for name, sel in (("converged_all", [t for t in table if t["converged"]]),
                      ("converged_redox_metal", [t for t in table if t["converged"] and t["redox"]])):
        dft = [t["dft"] for t in sel]
        out[name] = {"vs_label": stats(dft, [t["label"] for t in sel]),
                     "vs_regressor": stats(dft, [t["pred"] for t in sel]),
                     "members": [f"{t['jid']} {t['formula']}" for t in sel]}
    json.dump(out, open(HERE / "spotcheck_stats.json", "w"), indent=1)

    print("excluded (unconverged endpoint):", ", ".join(out["excluded_unconverged"]))
    for name in ("converged_all", "converged_redox_metal"):
        for ref in ("vs_label", "vs_regressor"):
            s = out[name][ref]
            print(f"{name:22s} {ref:13s} n={s['n']:2d}  MAE {s['MAE_V']:.2f}  ME {s['ME_V']:+.2f}  "
                  f"RMSE {s['RMSE_V']:.2f}  MAE(bias removed) {s['MAE_bias_removed_V']:.2f}  "
                  f"Spearman {s['spearman']:.2f}  window {s['window_agreement']}")
    print("wrote", HERE / "spotcheck_stats.json")


if __name__ == "__main__":
    main()
