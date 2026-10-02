"""Tier-1-as-front ranking after removing the label-set offset.

The tier-1 regressor is trained on the bmat_sign.py ALIGNN-FF labels, which sit
0.633 V below the Li_min.csv ALIGNN-FF screening values on average. This script
shifts the tier-1 predictions by that mean offset before applying the identical
screening filters, so that the window criterion is applied on the same voltage
scale as the tier-2 ranking, and reports the overlap.
"""
import json
import os
import sys

import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
import screen_cathode as sc  # noqa: E402

pred = pd.read_csv(os.path.join(HERE, "tier1_pool_predictions.csv"))
t2 = pd.read_csv(os.path.join(HERE, "..", "..", "cathode_candidates_ranked.csv"))
full = pd.read_csv(os.path.join(HERE, "tier1_front_ranking.csv"))
offset = (pred.v_tier2 - pred.v_tier1).mean()
print(f"mean offset tier2 - tier1 over pool: {offset:.3f} V")
# build pool frame with q_grav/ehull/formula from the uncorrected front file if present
cols = [c for c in full.columns if c in ("jid", "q_grav", "ehull", "formula", "max_voltage")]
qe = full[cols].drop_duplicates("jid") if "q_grav" in full.columns else None
if qe is None or len(qe) < len(pred):
    import pickle
    d = pickle.load(open("/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/jarvis_dft3d_cache.pkl", "rb"))
    df3 = pd.DataFrame(d)[["jid", "ehull", "formula", "atoms"]]
    df3 = df3[df3.jid.isin(pred.jid)].copy()
    df3["q_grav"] = df3.atoms.apply(sc.theoretical_grav_capacity)
    qe = df3[["jid", "ehull", "q_grav"]]
m = pred.merge(qe, on="jid", how="left")
m["avg_voltage"] = m.v_tier1 + offset
m["max_voltage"] = m.max_voltage_tier2
surv = sc.filter_cathode_candidates(m)
rk = sc.rank_candidates(surv)
rk.insert(0, "rank", range(1, len(rk) + 1))
common = set(rk.jid) & set(t2.jid)
a = rk.set_index("jid").loc[sorted(common), "rank"]
b = t2.set_index("jid").loc[sorted(common), "rank"]
res = {
    "offset_applied_V": round(offset, 3),
    "n_survivors_tier1_offset": len(rk),
    "n_survivors_tier2": len(t2),
    "n_common": len(common),
    "jaccard": len(common) / len(set(rk.jid) | set(t2.jid)),
    "spearman_common": spearmanr(a, b).correlation,
    "top12_overlap": len(set(rk.jid.head(12)) & set(t2.jid.head(12))),
    "top50_overlap": len(set(rk.jid.head(50)) & set(t2.jid.head(50))),
    "top12_tier1_offset": list(rk.jid.head(12)),
    "rank_LCO": int(rk[rk.jid == "JVASP-2017"]["rank"].iloc[0]) if (rk.jid == "JVASP-2017").any() else None,
    "rank_Li3FeO3": int(rk[rk.jid == "JVASP-117295"]["rank"].iloc[0]) if (rk.jid == "JVASP-117295").any() else None,
}
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(HERE, "front_ranking_offset.json"), "w"), indent=1)
rk.drop(columns=[c for c in rk.columns if c == "atoms"]).to_csv(os.path.join(HERE, "tier1_front_ranking_offset.csv"), index=False)
