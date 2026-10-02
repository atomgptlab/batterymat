"""Task 6 (R1.m4): tier-1 predictions for the five DFT benchmark cathodes."""
import json
import os
import sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict import ARCHIVE, MODEL_DIR, Tier1Predictor, load_cache  # noqa: E402

BENCH = {"JVASP-42723": "LFP LiFePO4", "JVASP-116897": "LMP LiMnPO4", "JVASP-141792": "LMO LiMn2O4",
         "JVASP-144791": "NMC Li4Mn3Co2Ni3O16", "JVASP-2017": "LCO LiCoO2"}
DFT_AVG = {"JVASP-42723": 3.60, "JVASP-116897": 3.91, "JVASP-141792": 4.08, "JVASP-144791": 4.40, "JVASP-2017": 4.18}
EXP = {"JVASP-42723": "3.45", "JVASP-116897": "~4.1", "JVASP-141792": "4.1", "JVASP-144791": "~3.7", "JVASP-2017": "3.9-4.2"}

ids = json.load(open(os.path.join(MODEL_DIR, "ids_train_val_test.json")))
part = {j: k.replace("id_", "") for k, v in ids.items() for j in v}
summ = pd.read_csv(os.path.join(ARCHIVE, "summary.csv"))
summ["avg_voltage_V"] = pd.to_numeric(summ["avg_voltage_V"], errors="coerce")
label = dict(zip(summ["jid"], summ["avg_voltage_V"]))
li = pd.read_csv(os.path.join(HERE, "..", "..", "..", "average_voltage", "Li_min.csv"))
li = li[li["name"].str.startswith("Li_")]
li_v = {n.split("_")[1]: v for n, v in zip(li["name"], li["avg_voltage"])}
# archived predictions: Test_results.json is per-structure with ids. Val_results
# .json is batched (32/batch, drop_last) in id_val order (verified below on the
# targets). Train_results.json comes from a shuffled loader, so archived
# per-id train predictions cannot be recovered; the column is NaN for train.
res = {}
for r in json.load(open(os.path.join(MODEL_DIR, "Test_results.json"))):
    res[r["id"]] = r["pred_out"]
idprop = {d["jid"]: d["target"] for d in json.load(open(os.path.join(ARCHIVE, "id_prop.json")))}
batches = json.load(open(os.path.join(MODEL_DIR, "Val_results.json")))
vt = [v for b in batches for v in b["target_out"]]; vp = [v for b in batches for v in b["pred_out"]]
assert max(abs(idprop[j] - t) for j, t in zip(ids["id_val"], vt)) < 1e-4
for j, v in zip(ids["id_val"], vp):
    res[j] = v
cache = load_cache()
p = Tier1Predictor()
rows = []
for j, name in BENCH.items():
    rows.append({"jid": j, "material": name, "formula": cache[j]["formula"], "n_atoms": len(cache[j]["atoms"]["elements"]),
                 "partition": part.get(j, "none"), "v_tier1_V": p.predict_one(cache[j]["atoms"]),
                 "v_tier1_archived_V": res.get(j), "v_label_summary_V": label.get(j), "v_tier2_Li_min_V": li_v.get(j),
                 "v_dft_avg_V": DFT_AVG[j], "v_exp_V": EXP[j]})
df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "benchmarks_tier1.csv"), index=False)
print(df.to_string(index=False))
