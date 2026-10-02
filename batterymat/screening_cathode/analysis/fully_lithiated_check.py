"""Proxy check of the fully-lithiated-input assumption (reviewer R2.M9).

For every entry of the lithium screening pool, look for another JARVIS-DFT
entry with the same host stoichiometry (identical non-Li composition after
normalisation) but a higher Li:host ratio. If one exists, the pool entry is
not the most lithiated form of its host known to the database.
"""
import pickle
import re
from collections import defaultdict
from fractions import Fraction

import pandas as pd

CACHE = "/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/jarvis_dft3d_cache.pkl"
POOL = "/Users/jaelee/Desktop/work/batterymat_jae/batterymat_jae/average_voltage/Li_min.csv"


def parse(f):
    d = defaultdict(int)
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        d[el] += int(n) if n else 1
    return d


def host_key(comp):
    host = {k: v for k, v in comp.items() if k != "Li"}
    if not host:
        return None, None
    from math import gcd
    from functools import reduce
    g = reduce(gcd, host.values())
    key = tuple(sorted((k, v // g) for k, v in host.items()))
    li_ratio = Fraction(comp.get("Li", 0), g)
    return key, li_ratio


d = pickle.load(open(CACHE, "rb"))
by_host = defaultdict(list)
for x in d:
    key, r = host_key(parse(x["formula"]))
    if key is not None and r > 0:
        by_host[key].append((r, x["jid"]))
li = pd.read_csv(POOL)
li = li[li.name.str.startswith("Li_")]
li["jid"] = li.name.str.split("_").str[1]
li["formula"] = li.name.str.split("_").str[2].str.replace(".json", "", regex=False)
n_viol = 0
viol = []
for _, r in li.iterrows():
    key, ratio = host_key(parse(r.formula))
    if key is None:
        continue
    higher = [j for rr, j in by_host[key] if rr > ratio]
    if higher:
        n_viol += 1
        viol.append((r.jid, r.formula, higher[:3]))
print(f"pool {len(li)}; entries with a more-lithiated same-host entry in JARVIS-DFT: {n_viol} ({100*n_viol/len(li):.1f}%)")
bench = {"JVASP-42723", "JVASP-116897", "JVASP-141792", "JVASP-144791", "JVASP-2017"}
print("benchmarks flagged:", [v for v in viol if v[0] in bench])
pd.DataFrame(viol, columns=["jid", "formula", "more_lithiated_jids"]).to_csv("fully_lithiated_violations.csv", index=False)
