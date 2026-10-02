"""Build composition-family-disjoint (grouped) train/test splits for the
ALIGNN average-voltage regressor (reviewer R1.M7 / R2.M7 / R2.m2).

Reads the archived training set (id_prop.json: list of {jid, atoms, target}),
assigns every structure a chemistry family, and writes k fold directories, each
holding id_prop.json (train, val, test concatenated in that order; the job sets
n_train/n_val/n_test and keep_data_order=true) plus the three parts separately. Two grouping schemes are produced:

  redox   : dominant redox-active transition metal (Ti V Cr Mn Fe Co Ni Cu Nb Mo
            Ru Rh W), 'none' otherwise      -> leave-one-family-out
  anion   : anion class (O-only oxide, F-containing, S/Se, P-O polyanion, other)
            -> leave-one-class-out

Nothing is trained here. Run on the cluster with job_grouped_splits.sh.

Usage:
  python make_grouped_splits.py --id-prop /path/to/Li_250/id_prop.json --out folds/
"""
import argparse
import json
import os
import random
import re
from collections import Counter, defaultdict

REDOX = ["Mn", "Fe", "Co", "Ni", "V", "Cr", "Ti", "Cu", "Nb", "Mo", "Ru", "Rh", "W"]


def elements_of(atoms):
    return atoms["elements"]


def redox_family(els):
    c = Counter(e for e in els if e in REDOX)
    return c.most_common(1)[0][0] if c else "none"


def anion_family(els):
    s = set(els)
    if "F" in s:
        return "fluoride"
    if "S" in s or "Se" in s:
        return "chalcogenide"
    if "P" in s and "O" in s:
        return "phosphate"
    if "O" in s:
        return "oxide"
    return "other"


def write_fold(out, name, train, val, test):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    for tag, rows in (("train", train), ("val", val), ("test", test)):
        json.dump(rows, open(os.path.join(d, f"id_prop_{tag}.json"), "w"))
    # train_alignn reads a single id_prop.json; with keep_data_order=true and
    # explicit n_train/n_val/n_test it takes the first n_train rows as train,
    # the next n_val as val and the last n_test as test.
    json.dump(train + val + test, open(os.path.join(d, "id_prop.json"), "w"))
    json.dump({"n_train": len(train), "n_val": len(val), "n_test": len(test),
               "test_families": sorted({r["family"] for r in test})},
              open(os.path.join(d, "fold_info.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id-prop", required=True)
    ap.add_argument("--out", default="folds")
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--val-frac", type=float, default=0.1)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    data = json.load(open(args.id_prop))
    for r in data:
        els = elements_of(r["atoms"])
        r["redox"] = redox_family(els)
        r["anion"] = anion_family(els)
    for scheme in ("redox", "anion"):
        groups = defaultdict(list)
        for r in data:
            groups[r[scheme]].append(r)
        print(f"[{scheme}] families:", {k: len(v) for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))})
        for fam, rows in groups.items():
            if len(rows) < 50:
                continue  # too small to be a meaningful held-out family
            test = [dict(r, family=fam) for r in rows]
            rest = [dict(r, family=r[scheme]) for r in data if r[scheme] != fam]
            rng.shuffle(rest)
            n_val = int(round(args.val_frac * len(rest)))
            write_fold(args.out, f"{scheme}_holdout_{fam}", rest[n_val:], rest[:n_val], test)
    # random k-fold (k=5) for MAE stability, same seed
    idx = list(range(len(data)))
    rng.shuffle(idx)
    k = 5
    for f in range(k):
        test = [dict(data[i], family="random") for i in idx[f::k]]
        rest = [dict(data[i], family="random") for j, i in enumerate(idx) if j % k != f]
        n_val = int(round(args.val_frac * len(rest)))
        write_fold(args.out, f"random_fold{f}", rest[n_val:], rest[:n_val], test)
    print("folds written to", args.out)


if __name__ == "__main__":
    main()
