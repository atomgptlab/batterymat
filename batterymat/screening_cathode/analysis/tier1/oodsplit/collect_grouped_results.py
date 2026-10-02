"""Collect per-fold held-out metrics from ALIGNN outputs (prediction_results_test_set.csv)."""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "folds"
rows = []
for d in sorted(glob.glob(os.path.join(root, "*"))):
    f = os.path.join(d, "out", "prediction_results_test_set.csv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f)
    y, p = df.iloc[:, 1].values, df.iloc[:, 2].values
    info = json.load(open(os.path.join(d, "fold_info.json")))
    mae = np.abs(y - p).mean()
    r2 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rows.append({"fold": os.path.basename(d), "n_test": len(y), "test_families": ";".join(info["test_families"]),
                 "MAE_V": round(mae, 4), "R2": round(r2, 4), "MAD_V": round(np.abs(y - y.mean()).mean(), 4)})
out = pd.DataFrame(rows)
print(out.to_csv(index=False))
