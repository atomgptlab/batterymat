"""Task 2: tier-1 predictions over the full Li_ screening pool (Li_min.csv).

Writes tier1_pool_predictions.csv with columns
  jid, formula, v_tier1, v_tier2 (Li_min.csv avg_voltage), v_label
  (summary.csv avg_voltage_V, NaN if absent), partition (train/val/test/none).
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict import ARCHIVE, MODEL_DIR, Tier1Predictor, load_cache  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LI_MIN = os.path.join(REPO, "average_voltage", "Li_min.csv")
OUT_CSV = os.path.join(HERE, "tier1_pool_predictions.csv")


def pool_frame():
    li = pd.read_csv(LI_MIN)
    li = li[li["name"].str.startswith("Li_")].copy()
    li["jid"] = li["name"].str.split("_").str[1]
    li = li.drop_duplicates("jid").reset_index(drop=True)
    return li[["jid", "avg_voltage", "max_voltage"]].rename(columns={"avg_voltage": "v_tier2", "max_voltage": "max_voltage_tier2"})


def main(limit=None):
    li = pool_frame()
    print(f"Li_ pool: {len(li)} unique JIDs")
    cache = load_cache()
    summ = pd.read_csv(os.path.join(ARCHIVE, "summary.csv"))
    summ["avg_voltage_V"] = pd.to_numeric(summ["avg_voltage_V"], errors="coerce")
    label = dict(zip(summ["jid"], summ["avg_voltage_V"]))
    ids = json.load(open(os.path.join(MODEL_DIR, "ids_train_val_test.json")))
    part = {}
    for k, name in (("id_train", "train"), ("id_val", "val"), ("id_test", "test")):
        for j in ids[k]:
            part[j] = name
    if limit:
        li = li.head(limit)
    missing = [j for j in li["jid"] if j not in cache]
    print(f"JIDs missing from JARVIS cache: {len(missing)} {missing[:10]}")
    p = Tier1Predictor()
    t0 = time.time()
    preds, formulas = [], []
    for i, j in enumerate(li["jid"]):
        if j in cache:
            try:
                preds.append(p.predict_one(cache[j]["atoms"]))
            except Exception as exc:
                print(f"[warn] {j} failed: {exc}")
                preds.append(np.nan)
            formulas.append(cache[j]["formula"])
        else:
            preds.append(np.nan)
            formulas.append("")
        if (i + 1) % 100 == 0:
            el = time.time() - t0
            print(f"{i+1}/{len(li)} {el:.0f}s ({el/(i+1):.3f} s/struct, ETA {el/(i+1)*(len(li)-i-1)/60:.1f} min)", flush=True)
    li["formula"] = formulas
    li["v_tier1"] = preds
    li["v_label"] = li["jid"].map(label)
    li["partition"] = li["jid"].map(part).fillna("none")
    li = li[["jid", "formula", "v_tier1", "v_tier2", "v_label", "partition", "max_voltage_tier2"]]
    li.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(li)} rows) in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
