#!/usr/bin/env python3
"""BatteryMat TODO item 2: ranking/enrichment analysis from existing predictions.

Reproduces the paper's cathode-window pool (expected n=2,890), ranks it by the
tier-2 (ALIGNN-FF) screening quantities, and evaluates ranking against the five
DFT-validated cathodes. Optionally merges tier-1 single-shot predictions
(--tier1-csv, columns: jid,pred_v) once pulled from the cluster, and then
reports Spearman(tier-1, tier-2) over the full pool.

Run:  python3 02_enrichment.py
Output: out/enrichment_tables.md
"""

import argparse
import csv
import json
import os
import re

import numpy as np

REPO = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
)
LI_MIN_CSV = os.path.join(REPO, "average_voltage/Li_min.csv")
RANKED_CSV = os.path.join(REPO, "screening_cathode/cathode_candidates_ranked.csv")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# DFT-validated benchmarks: JVASP id -> (label, DFT avg V, exp avg V)
# DFT/exp values: fill/verify from the paper's Table 1 before pasting results.
DFT_ANCHORS = {
    "JVASP-42723": ("LFP", None, 3.45),
    "JVASP-116897": ("LMP", None, 4.1),
    "JVASP-141792": ("LMO", None, 4.15),
    "JVASP-144791": ("NMC", None, 3.8),
    "JVASP-2017": ("LCO", None, 3.9),
}

PAPER_POOL_N = 2890
V_MAX, V_MIN = 5.5, 1.0


def load_pool():
    rows, seen = [], set()
    with open(LI_MIN_CSV) as f:
        for r in csv.DictReader(f):
            m = re.search(r"(JVASP-\d+)", r.get("name", ""))
            if not m:
                continue
            try:
                avg_v = float(r["avg_voltage"])
                max_v = float(r["max_voltage"])
                profile = [float(x) for x in r["voltage_profile"].split(";") if x]
                min_v = min(profile) if profile else avg_v
            except (KeyError, ValueError):
                continue
            jid = m.group(1)
            if jid in seen:
                continue
            seen.add(jid)
            rows.append({"jid": jid, "avg_v": avg_v, "max_v": max_v,
                         "min_v": min_v, "name": r.get("name", "")})
    return rows


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    denom = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float((ra * rb).sum() / denom) if denom else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier1-csv", default=None,
                    help="optional csv with columns jid,pred_v (tier-1 surrogate predictions)")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    lines = ["# Ranking / enrichment analysis (existing predictions only)", ""]

    pool = load_pool()
    win = [r for r in pool if r["max_v"] <= V_MAX and r["min_v"] > V_MIN]
    check = ("OK" if len(win) == PAPER_POOL_N else
             f"paper reports {PAPER_POOL_N}; delta likely an additional completion "
             f"constraint or pool-snapshot difference - RECONCILE BEFORE QUOTING")
    lines.append(f"Cathode-window pool: **{len(win)}** unique materials "
                 f"(max V <= {V_MAX}, min V > {V_MIN}, dedup by jid)  [{check}]")
    lines.append("")

    # rank pool by tier-2 average voltage (descending) and locate the DFT anchors
    ranked = sorted(win, key=lambda r: -r["avg_v"])
    jid_rank = {r["jid"]: i + 1 for i, r in enumerate(ranked)}
    n = len(ranked)
    lines.append("## Where the five DFT-validated cathodes fall in the tier-2 ranking")
    lines.append("| Cathode | Tier-2 rank (of pool) | Percentile | Tier-2 avg V | Exp avg V |")
    lines.append("|---|---|---|---|---|")
    anchor_rows = []
    for jid, (label, dft_v, exp_v) in DFT_ANCHORS.items():
        rk = jid_rank.get(jid)
        r = next((x for x in ranked if x["jid"] == jid), None)
        if rk is None:
            lines.append(f"| {label} ({jid}) | not in window pool | - | - | {exp_v} |")
            continue
        anchor_rows.append((label, r["avg_v"], exp_v))
        lines.append(f"| {label} ({jid}) | {rk} | {100*(1-rk/n):.1f} | {r['avg_v']:.2f} | {exp_v} |")
    lines.append("")

    if len(anchor_rows) >= 3:
        rho = spearman([a[1] for a in anchor_rows], [a[2] for a in anchor_rows])
        lines.append(f"Spearman(tier-2 avg V, experimental avg V) over the "
                     f"{len(anchor_rows)} in-pool anchors: **{rho:.2f}** "
                     f"(N is small - state this honestly; it motivates the prospective run).")
        lines.append("")

    if args.tier1_csv:
        t1 = {}
        with open(args.tier1_csv) as f:
            for r in csv.DictReader(f):
                t1[r["jid"]] = float(r["pred_v"])
        both = [(t1[r["jid"]], r["avg_v"]) for r in ranked if r["jid"] in t1]
        rho = spearman([b[0] for b in both], [b[1] for b in both])
        lines.append(f"## Tier-1 vs tier-2 ranking agreement")
        lines.append(f"Spearman over {len(both)} pool structures with tier-1 predictions: **{rho:.3f}**")
        lines.append("")
    else:
        lines.append("_Tier-1 single-shot predictions not merged (no --tier1-csv). Pull the "
                     "surrogate's pool predictions from the training cluster to add "
                     "Spearman(tier-1, tier-2) over the full pool._")
        lines.append("")

    out = os.path.join(OUT_DIR, "enrichment_tables.md")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
