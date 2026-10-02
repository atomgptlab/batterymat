"""Campaign C: DFT endpoint spot-check of the tier-1 ALIGNN voltage labels.

Selects 30 structures from the tier-1 held-out test partition that are also
in the Li screening pool, stratified over the label voltage range, and writes
two primitive-cell VASP relaxations per structure (fully lithiated and fully
delithiated) under ``dft_inputs/spotcheck/<JID>/{lithiated,delithiated}/``.

Selection filters (in order):
  1. jid in tier-1 held-out test split (ids_train_val_test.json, key ``id_test``)
  2. jid in the Li screening pool (average_voltage/Li_min.csv rows named ``Li_*``)
  3. primitive cell <= 20 atoms (JARVIS dft_3d snapshot)
  4. ehull <= 0.10 eV/atom
  5. no f-block elements
  6. stratified sample: 6 bins of 1 V over 0-6 V of the label voltage
     (summary.csv ``avg_voltage_V``), 5 per bin, fixed seed; a short bin is
     filled from the nearest neighbouring bins.

Average DFT voltage from the two endpoints (same convention as
``dft_prep.compute_voltage_curve``):

    V = (E_delithiated - E_lithiated) / n_Li + E_Li_metal

Inputs reuse the dft_prep.py INCAR template (ISIF=3 for both endpoints),
auto-selected functional (optB88-vdW for layered spacegroups, PBE otherwise),
DFT+U and AFM MAGMOM exactly as dft_prep generates them. KPOINTS are
Gamma-centred with n_i = max(1, ceil(21 A / |a_i|)), which reproduces the
3x3x3 mesh dft_prep uses at its 7 A supercell threshold.

Usage (from screening_cathode/analysis/):
    python spotcheck_prepare.py [--archive DIR] [--db-cache PKL] [--seed N] [--n 30]
"""
import argparse
import json
import math
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SC_DIR = HERE.parent                      # screening_cathode/
REPO = SC_DIR.parents[1]                  # repository root
sys.path.insert(0, str(SC_DIR))
import dft_prep as dp  # noqa: E402

from jarvis.core.atoms import Atoms  # noqa: E402
from jarvis.io.vasp.inputs import Poscar  # noqa: E402

DEFAULT_ARCHIVE = Path("/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250")
F_BLOCK = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
              "Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split())
KPOINT_LENGTH = 21.0   # A; 21/7 = 3 -> matches dft_prep's 3x3x3 at the 7 A threshold


def _jid_of(name: str):
    m = re.search(r"(JVASP-\d+)", str(name))
    return m.group(1) if m else None


def kmesh_for(atoms):
    norms = np.linalg.norm(np.array(atoms.lattice_mat), axis=1)
    return " ".join(str(max(1, math.ceil(KPOINT_LENGTH / n))) for n in norms)


def delithiate(atoms):
    keep = [i for i, el in enumerate(atoms.elements) if el != "Li"]
    return Atoms(
        lattice_mat=atoms.lattice_mat,
        elements=[atoms.elements[i] for i in keep],
        coords=[atoms.frac_coords[i] for i in keep],
        cartesian=False,
    )


def write_endpoint(directory: Path, atoms, functional: str) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    elements = list(dict.fromkeys(atoms.elements))
    n_atoms = len(atoms.elements)
    nsw = min(300, max(200, n_atoms * 2))
    magmom = dp._magmom_string(dp._afm_magmom(atoms))
    Poscar(atoms).write_file(str(directory / "POSCAR"))
    dp.write_relax_incar(directory, elements=elements, isif=3, nsw=nsw,
                         magmom_str=magmom, functional=functional)
    mesh = kmesh_for(atoms)
    dp.write_kpoints(directory, mesh=mesh)
    dp.write_potcar_spec(directory, elements)
    return {"n_atoms": n_atoms, "kmesh": mesh, "nsw": nsw}


def select(args):
    archive = Path(args.archive)
    ids_repo = REPO / "batterymat_jae" / "average_voltage" / "Li_250_ids_train_val_test.json"
    ids_path = ids_repo if ids_repo.exists() else archive / "voltage" / "ids_train_val_test.json"
    ids = json.loads(Path(ids_path).read_text())
    test_ids = set(ids["id_test"])
    print(f"tier-1 test split: {len(test_ids)} ids from {ids_path}")

    li_min = pd.read_csv(REPO / "batterymat_jae" / "average_voltage" / "Li_min.csv")
    li_min = li_min[li_min["name"].str.startswith("Li_")].copy()
    li_min["jid"] = li_min["name"].map(_jid_of)
    li_pool = li_min.set_index("jid")
    print(f"Li screening pool: {len(li_pool)} rows")

    labels = pd.read_csv(archive / "summary.csv").set_index("jid")
    print(f"labels (summary.csv avg_voltage_V): {len(labels)} rows")

    rows = pickle.load(open(args.db_cache, "rb"))
    db = {r["jid"]: r for r in rows}

    eligible = []
    for jid in sorted(test_ids & set(li_pool.index)):
        r = db.get(jid)
        if r is None or jid not in labels.index:
            continue
        els = r["atoms"]["elements"]
        ehull = r["ehull"]
        if len(els) > 20 or not isinstance(ehull, (int, float)) or ehull > 0.10:
            continue
        if set(els) & F_BLOCK:
            continue
        v = float(labels.loc[jid, "avg_voltage_V"])
        if not (0.0 <= v <= 6.0):
            continue
        spg = int(r["spg_number"])
        eligible.append({
            "jid": jid, "formula": r["formula"], "n_atoms_prim": len(els),
            "n_li": sum(1 for e in els if e == "Li"),
            "ehull_eV_atom": float(ehull), "spg_number": spg,
            "label_voltage_V": v,
            "li_min_avg_voltage_V": float(li_pool.loc[jid, "avg_voltage"]),
            "bin": min(int(v // 1.0), 5),
        })
    elig = pd.DataFrame(eligible)
    print(f"eligible after filters: {len(elig)}; per-bin counts: "
          f"{elig['bin'].value_counts().sort_index().to_dict()}")

    rng = np.random.default_rng(args.seed)
    per_bin = args.n // 6
    chosen = []
    remaining = {b: list(rng.permutation(elig.index[elig["bin"] == b])) for b in range(6)}
    for b in range(6):
        take = remaining[b][:per_bin]
        remaining[b] = remaining[b][per_bin:]
        chosen.extend((i, b, "primary") for i in take)
        short = per_bin - len(take)
        # fill a short bin from the nearest neighbouring bins
        for dist in range(1, 6):
            if short == 0:
                break
            for nb in (b - dist, b + dist):
                if 0 <= nb < 6 and short > 0 and remaining[nb]:
                    i = remaining[nb].pop(0)
                    chosen.append((i, b, f"fill_from_bin_{nb}"))
                    short -= 1
    sel = elig.loc[[i for i, _, _ in chosen]].copy()
    sel["target_bin"] = [b for _, b, _ in chosen]
    sel["selection"] = [s for _, _, s in chosen]
    sel = sel.sort_values(["target_bin", "label_voltage_V"]).reset_index(drop=True)
    return sel, db


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--archive", default=str(DEFAULT_ARCHIVE),
                    help="Li_250 archive dir (summary.csv, voltage/ids_train_val_test.json)")
    ap.add_argument("--db-cache", default=str(DEFAULT_ARCHIVE / "jarvis_dft3d_cache.pkl"),
                    help="JARVIS dft_3d snapshot pickle")
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--out", default=str(SC_DIR / "dft_inputs" / "spotcheck"))
    args = ap.parse_args()

    sel, db = select(args)
    out_root = Path(args.out)
    manifest = []
    for _, row in sel.iterrows():
        jid = row["jid"]
        r = db[jid]
        atoms = Atoms.from_dict(r["atoms"])
        layered = int(r["spg_number"]) in dp._LAYERED_SPACEGROUPS
        functional = "optb88vdw" if layered else "pbe"
        e_li = dp._E_LI_METAL[functional]
        lith = write_endpoint(out_root / jid / "lithiated", atoms, functional)
        delith = write_endpoint(out_root / jid / "delithiated", delithiate(atoms), functional)
        manifest.append({
            "jid": jid, "formula": r["formula"], "spg_number": int(r["spg_number"]),
            "ehull_eV_atom": float(r["ehull"]), "n_li": int(row["n_li"]),
            "n_atoms_lithiated": lith["n_atoms"], "n_atoms_delithiated": delith["n_atoms"],
            "functional": functional, "layered": layered, "e_li_metal": e_li,
            "isif": 3, "kmesh": lith["kmesh"],
            "label_voltage_V": float(row["label_voltage_V"]),
            "voltage_formula": "V = (E_delithiated - E_lithiated) / n_li + e_li_metal",
        })
        print(f"{jid:13s} {r['formula']:18s} bin={row['target_bin']} label={row['label_voltage_V']:.2f} V "
              f"{functional:9s} nat={lith['n_atoms']}/{delith['n_atoms']} k={lith['kmesh']}")

    sel["functional"] = [m["functional"] for m in manifest]
    sel["e_li_metal"] = [m["e_li_metal"] for m in manifest]
    sel["rationale"] = (
        "tier-1 held-out test id; in Li screening pool; primitive cell <=20 atoms; "
        "ehull<=0.10 eV/atom; no f-block; stratified 5 per 1 V bin over 0-6 V "
        f"(seed {args.seed})"
    )
    sel_path = HERE / "spotcheck_selection.csv"
    sel.to_csv(sel_path, index=False)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "spotcheck_manifest.json").write_text(json.dumps({
        "campaign": "C: DFT endpoint spot-check of tier-1 ALIGNN voltage labels",
        "seed": args.seed, "n": len(manifest),
        "label_source": "Li_250/summary.csv avg_voltage_V",
        "voltage_formula": "V = (E_delithiated - E_lithiated) / n_li + e_li_metal",
        "e_li_metal": dp._E_LI_METAL,
        "structures": manifest,
    }, indent=2) + "\n")
    print(f"wrote {sel_path} and {out_root / 'spotcheck_manifest.json'}")


if __name__ == "__main__":
    main()
