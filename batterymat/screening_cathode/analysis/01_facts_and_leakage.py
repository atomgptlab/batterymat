#!/usr/bin/env python3
"""BatteryMat TODO items 1+4: dataset facts and the five-cathode leakage check.

Outputs out/facts_and_leakage.md answering:
  - training-pool size (label-generation set from screening_anode/Li.csv) vs the paper's 7,610
  - whether each of the five benchmark cathodes is IN the label pool (SP5Z §1a)
  - structural near-duplicates of the five cathodes inside the pool
    (reduced-formula bucket -> pymatgen StructureMatcher, vs jarvis dft_3d structures)

Run:  python3 01_facts_and_leakage.py
"""

import csv
import json
import os
import re
import sys

REPO = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
)
LI_CSV = os.path.join(REPO, "screening_anode/Li.csv")
LI_MIN_CSV = os.path.join(REPO, "average_voltage/Li_min.csv")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
# jarvis-tools corpus downloads land here on first run (override via env)
JARVIS_CACHE_DIR = os.environ.get(
    "JARVIS_CACHE_DIR", os.path.join(OUT_DIR, "jarvis_cache")
)

FIVE = {
    "JVASP-42723": "LFP (LiFePO4)",
    "JVASP-116897": "LMP (LiMnPO4)",
    "JVASP-141792": "LMO (LiMn2O4)",
    "JVASP-144791": "NMC (Li[Ni,Mn,Co]O2)",
    "JVASP-2017": "LCO (LiCoO2)",
}


def jids_from_csv(path):
    jids = set()
    with open(path) as f:
        for row in csv.DictReader(f):
            jid = row.get("jid", "")
            if not re.match(r"JVASP-\d+$", jid or ""):
                m = re.search(r"(JVASP-\d+)", row.get("name", ""))
                jid = m.group(1) if m else None
            if jid:
                jids.add(jid)
    return jids


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pool = jids_from_csv(LI_CSV)
    pool_min = jids_from_csv(LI_MIN_CSV)
    lines = ["# Dataset facts + five-cathode leakage check", ""]
    lines.append(f"- Label pool (screening_anode/Li.csv): **{len(pool)}** unique JVASP ids "
                 f"(paper states 7,610 structures in the training pool - reconcile the delta: "
                 f"attempted vs completed label generation).")
    lines.append(f"- Filtered pool (average_voltage/Li_min.csv): {len(pool_min)} ids.")
    lines.append("")

    # membership check
    rows = [(jid, name, "YES - in label pool" if jid in pool else "no")
            for jid, name in FIVE.items()]
    lines.append("## Five benchmark cathodes vs the label pool (exact membership)")
    lines.append("| JVASP id | Material | In training/label pool? |")
    lines.append("|---|---|---|")
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} |")
    lines.append("")

    # near-duplicate check against dft_3d structures sharing a reduced formula
    print("fetching dft_3d for structural near-dup check (cached if already downloaded)...")
    cwd = os.getcwd()
    os.makedirs(JARVIS_CACHE_DIR, exist_ok=True)
    os.chdir(JARVIS_CACHE_DIR)
    try:
        from jarvis.db.figshare import data as jdata
        entries = jdata("dft_3d")
    finally:
        os.chdir(cwd)
    by_jid = {e.get("jid"): e for e in entries}

    from pymatgen.core import Composition
    from jarvis.core.atoms import Atoms

    try:
        from pymatgen.core.structure_matcher import StructureMatcher
    except ImportError:
        from pymatgen.analysis.structure_matcher import StructureMatcher
    near_m = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True)

    def rf(e):
        try:
            return Composition(e["formula"]).reduced_formula
        except Exception:
            return None

    lines.append("## Structural near-duplicates of the five cathodes inside the label pool")
    lines.append("| Benchmark cathode | Same reduced formula in pool (n) | Near-duplicates (StructureMatcher) |")
    lines.append("|---|---|---|")
    pool_entries = [by_jid[j] for j in pool if j in by_jid]
    for jid, name in FIVE.items():
        e = by_jid.get(jid)
        if e is None:
            lines.append(f"| {name} | not found in dft_3d | - |")
            continue
        target_rf = rf(e)
        cands = [c for c in pool_entries if rf(c) == target_rf and c["jid"] != jid]
        sb = Atoms.from_dict(e["atoms"]).pymatgen_converter()
        dups = []
        for c in cands:
            try:
                if near_m.fit(sb, Atoms.from_dict(c["atoms"]).pymatgen_converter()):
                    dups.append(c["jid"])
            except Exception:
                continue
        lines.append(f"| {name} ({jid}) | {len(cands)} | {', '.join(dups) if dups else 'none'} |")
    lines.append("")
    lines.append(
        "Interpretation for the rebuttal: exact membership answers whether the five cathodes "
        "were in-sample for the *label pool*; whether they were in the surrogate's train vs "
        "held-out TEST split additionally requires the split record from the training run "
        "(on the cluster - see README)."
    )

    out = os.path.join(OUT_DIR, "facts_and_leakage.md")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
