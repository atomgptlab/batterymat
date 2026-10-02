"""Task 3: tier-1-as-front ranking vs the tier-2 (ALIGNN-FF) ranking.

Identical filters/score to screen_cathode.py (imported live, so any filter
added there, e.g. the redox-metal requirement, applies here too), with the tier-1 predicted
average voltage replacing avg_voltage. q_grav and ehull are taken from
cathode_candidates_ranked.csv where the JID is present, otherwise recomputed
from the archived JARVIS cache (theoretical_grav_capacity + ehull). The
max_voltage filter keeps the tier-2 (ALIGNN-FF) max_voltage.
Outputs: tier1_front_ranking.csv, front_ranking_comparison.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
from predict import load_cache  # noqa: E402
import screen_cathode as sc  # noqa: E402

RANKED = os.path.abspath(os.path.join(HERE, "..", "..", "cathode_candidates_ranked.csv"))
LCO = "JVASP-2017"
TOPS = ["JVASP-117419", "JVASP-141543", "JVASP-117295", "JVASP-116849", "JVASP-154749", "JVASP-96563", "JVASP-95533", "JVASP-97428"]
CONTROLS = ["JVASP-80802", "JVASP-77457", "JVASP-85076", "JVASP-81829"]


def failed_filters(row, vcol):
    f = []
    v = row[vcol]
    if not (sc.FILTERS["avg_voltage_min"] <= v <= sc.FILTERS["avg_voltage_max"]):
        f.append(f"{vcol}={v:.2f}")
    if not row["q_grav"] > sc.FILTERS["q_grav_min"]:
        f.append(f"q_grav={row['q_grav']:.1f}")
    if not row["ehull"] <= sc.FILTERS["ehull_max"]:
        f.append(f"ehull={row['ehull']:.3f}")
    if not row["max_voltage"] <= sc.FILTERS["max_voltage_max"]:
        f.append(f"max_voltage={row['max_voltage']:.2f}")
    if sc.FILTERS.get("require_redox_metal") and not sc.has_redox_metal(str(row["formula"])):
        f.append("no_redox_metal")
    return ";".join(f) if f else "pass"


def main():
    pool = pd.read_csv(os.path.join(HERE, "tier1_pool_predictions.csv"))
    ranked = pd.read_csv(RANKED)
    cache = load_cache()
    # q_grav / ehull: from the tier-2 ranking when present, else from cache
    rq = dict(zip(ranked["jid"], ranked["q_grav"])); re_ = dict(zip(ranked["jid"], ranked["ehull"]))
    q, e, src = [], [], []
    for j in pool["jid"]:
        if j in rq:
            q.append(rq[j]); e.append(re_[j]); src.append("ranked_csv")
        else:
            ent = cache[j]
            q.append(sc.theoretical_grav_capacity(ent["atoms"]))
            e.append(pd.to_numeric(ent["ehull"], errors="coerce")); src.append("cache")
    pool["q_grav"] = q; pool["ehull"] = e; pool["qe_source"] = src
    pool = pool.rename(columns={"max_voltage_tier2": "max_voltage"})
    # consistency check of recomputed q_grav vs ranked.csv on the overlap
    chk = pool[pool["jid"].isin(rq)].copy()
    chk["q_cache"] = [sc.theoretical_grav_capacity(cache[j]["atoms"]) for j in chk["jid"]]
    chk["e_cache"] = [pd.to_numeric(cache[j]["ehull"], errors="coerce") for j in chk["jid"]]
    qdiff = float(np.nanmax(np.abs(chk["q_cache"] - chk["q_grav"]))); ediff = float(np.nanmax(np.abs(chk["e_cache"] - chk["ehull"])))

    # tier-2 ranking rebuilt from the same pool frame (must equal ranked.csv)
    t2 = pool.rename(columns={"v_tier2": "avg_voltage"})
    t2 = sc.rank_candidates(sc.filter_cathode_candidates(t2.dropna(subset=["avg_voltage", "ehull"])))
    t2.insert(0, "rank", range(1, len(t2) + 1))
    same = (len(t2) == len(ranked)) and list(t2["jid"]) == list(ranked["jid"])
    # tier-1 ranking
    t1 = pool.rename(columns={"v_tier1": "avg_voltage"}).dropna(subset=["avg_voltage", "ehull"])
    t1 = sc.rank_candidates(sc.filter_cathode_candidates(t1))
    t1.insert(0, "rank", range(1, len(t1) + 1))
    t1 = t1.rename(columns={"avg_voltage": "v_tier1"})
    t2rank = dict(zip(ranked["jid"], ranked["rank"]))
    t1["rank_tier2"] = t1["jid"].map(t2rank)
    t1 = t1[["rank", "jid", "formula", "v_tier1", "v_tier2", "v_label", "max_voltage", "q_grav", "ehull", "score",
             "partition", "rank_tier2", "qe_source"]]
    t1.to_csv(os.path.join(HERE, "tier1_front_ranking.csv"), index=False)

    s1, s2 = list(t1["jid"]), list(ranked["jid"])
    common = [j for j in s1 if j in t2rank]
    r1 = dict(zip(t1["jid"], t1["rank"]))
    rho = float(spearmanr([r1[j] for j in common], [t2rank[j] for j in common]).correlation) if len(common) > 2 else float("nan")
    # rank correlation over the union with non-survivors counted as worst rank? report only common.
    def overlap(k):
        a, b = set(s1[:k]), set(s2[:k]); return {"k": k, "n_common": len(a & b), "jaccard": len(a & b) / len(a | b)}

    # where the named JIDs land
    pool_i = pool.set_index("jid")
    def where(j):
        if j not in pool_i.index:
            return {"jid": j, "in_pool": False}
        row = pool_i.loc[j]
        return {"jid": j, "in_pool": True, "formula": row["formula"], "partition": row["partition"],
                "v_tier1": float(row["v_tier1"]), "v_tier2": float(row["v_tier2"]), "q_grav": float(row["q_grav"]),
                "ehull": float(row["ehull"]), "max_voltage": float(row["max_voltage"]),
                "rank_tier1": int(r1[j]) if j in r1 else None, "rank_tier2": int(t2rank[j]) if j in t2rank else None,
                "tier1_filters": failed_filters(row, "v_tier1"), "tier2_filters": failed_filters(row, "v_tier2")}
    out = {
        "pool_size": int(len(pool)), "n_pool_missing_v_tier1": int(pool["v_tier1"].isna().sum()),
        "n_qe_from_ranked_csv": int((pool["qe_source"] == "ranked_csv").sum()), "n_qe_from_cache": int((pool["qe_source"] == "cache").sum()),
        "max_abs_diff_q_grav_cache_vs_ranked": qdiff, "max_abs_diff_ehull_cache_vs_ranked": ediff,
        "tier2_ranking_reproduced_from_pool": bool(same), "n_tier2_ranked_csv": int(len(ranked)), "n_tier2_rebuilt": int(len(t2)),
        "n_survivors_tier1": int(len(t1)), "n_survivors_tier2": int(len(ranked)),
        "n_common_survivors": len(common), "jaccard_full": len(set(s1) & set(s2)) / len(set(s1) | set(s2)),
        "n_tier1_only": int(len(set(s1) - set(s2))), "n_tier2_only": int(len(set(s2) - set(s1))),
        "spearman_rank_common": rho,
        "overlap_top12": overlap(12), "overlap_top50": overlap(50),
        "top12_tier1": s1[:12], "top12_tier2": s2[:12],
        "survivors_by_partition_tier1": t1["partition"].value_counts().to_dict(),
        "LCO": where(LCO), "tops": [where(j) for j in TOPS], "controls": [where(j) for j in CONTROLS],
    }
    json.dump(out, open(os.path.join(HERE, "front_ranking_comparison.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items() if k not in ("tops", "controls", "LCO")}, indent=1))
    print("LCO", out["LCO"])
    for r in out["tops"] + out["controls"]:
        print(r)


if __name__ == "__main__":
    main()
