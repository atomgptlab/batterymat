#!/usr/bin/env python
"""Exhaustive Li-vacancy enumeration vs. the greedy ALIGNN-FF ranking of dft_prep.py.

Reviewer R1.M6: does the greedy heuristic in ``dft_prep.get_next_vacancy`` (at
each step score every remaining Li site by the ALIGNN-FF single-point energy of
the structure with that Li removed, remove the lowest, continue) find the true
minimum-energy vacancy configuration at every composition?  Tested by brute
force on two 8-Li cells: every one of the 2^8 = 256 subsets of Li sites is
evaluated with ALIGNN-FF.

Two levels:
  (a) single-point energies on the unrelaxed vacancy structures (what
      ``get_next_vacancy`` ranks on);
  (b) ALIGNN-FF relaxation (ions + cell: ASE FIRE on an ExpCellFilter, the
      same optimizer/filter combination ``alignn.ff.ff.ForceField.optimize_atoms``
      uses; fmax and the step cap are CLI arguments).  ``dft_prep.py`` never
      relaxes with ALIGNN-FF (relaxation is done by VASP), so level (b) is the
      ML analogue of the DFT relax.

Three delithiation paths are compared per composition n_Li = 8..0:
  * greedy_sp      -- greedy on single-point energies of the *unrelaxed* lattice
                      (exact reproduction of get_next_vacancy applied to the
                      rigid step_00 cell, including its first-minimum tie-break);
  * greedy_contcar -- faithful reproduction of the production loop: at each step
                      the candidates are single-point energies computed on the
                      *relaxed* structure of the previous step (the DFT CONTCAR
                      in production, the ALIGNN-FF relaxed structure here), the
                      winner is relaxed, and the loop continues;
  * global_min     -- per-composition minimum over all C(8, n) configurations,
                      at the single-point level and (if relaxed) at the relaxed
                      level.

Voltages use the ALIGNN-FF Li reference computed here with the same calculator
(bcc Li, JVASP-913 geometry: single-point value for level (a), FIRE-relaxed value
for level (b)).  For comparison the JSON also records ``unary_energy("Li")`` =
-0.925 eV/atom, which is the JARVIS OptB88-vdW DFT value screen.py uses, and the
DFT Li_sv references used by dft_prep.py.

Usage
-----
  python enumerate.py JVASP-2017   [--relax auto|all|minimal|none] [--workers 4]
  python enumerate.py JVASP-96563  --ref-poscar /path/to/POSCAR
  python enumerate.py --summary    # rebuild enumeration_summary.md + PNGs from JSON

Environment: ALIGNN-FF on this Mac needs KMP_DUPLICATE_LIB_OK=TRUE and a single
OpenMP thread per process (multithreaded torch+dgl segfaults); the script sets
both and parallelizes over spawned single-thread worker processes instead.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import pickle
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "1")
warnings.filterwarnings("ignore")

import numpy as np

HERE = Path(__file__).resolve().parent
SCREENING_CATHODE = HERE.parents[1]          # .../screening_cathode
sys.path.insert(0, str(SCREENING_CATHODE))

DEFAULT_CACHE = ("/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/"
                 "Li_250/jarvis_dft3d_cache.pkl")
E_LI_UNARY_SCREEN_PY = -0.925      # jarvis unary_energy("Li"), used by screen.py
E_LI_DFT = {"pbe": -1.9031, "optb88vdw": -0.9646}   # dft_prep._E_LI_METAL
GAP_WINDOW_MEV = 25.0
DEGENERACY_TOL_EV = 1e-4           # energies equal to within this are one "level"

# ----------------------------------------------------------------------------
# Structure helpers
# ----------------------------------------------------------------------------

def load_from_cache(jid: str, cache: str):
    from jarvis.core.atoms import Atoms
    rows = pickle.load(open(cache, "rb"))
    for r in rows:
        if r["jid"] == jid:
            return Atoms.from_dict(r["atoms"]), r
    raise KeyError(f"{jid} not in {cache}")


def build_supercell(jid: str, cache: str):
    """Reproduce exactly what dft_prep.generate_sequential_init builds."""
    import dft_prep
    prim, row = load_from_cache(jid, cache)
    dim = dft_prep._min_supercell_dim(prim)          # same helper, same default 7 A
    sup = prim if dim == [1, 1, 1] else prim.make_supercell(dim)
    return sup, dim, row


def find_reference_poscar(jid: str):
    for d in sorted((SCREENING_CATHODE / "dft_inputs").glob(f"{jid}*")):
        for p in sorted(d.glob("supercell_*/step_00_*/POSCAR")):
            return p
    return None


def compare_to_poscar(atoms, poscar_path):
    from jarvis.io.vasp.inputs import Poscar
    ref = Poscar.from_file(str(poscar_path)).atoms
    if ref.elements != atoms.elements:
        return {"path": str(poscar_path), "match": False, "reason": "element list differs"}
    dl = float(np.abs(np.array(ref.lattice_mat) - np.array(atoms.lattice_mat)).max())
    df = np.array(ref.frac_coords) - np.array(atoms.frac_coords)
    df = float(np.abs(df - np.round(df)).max())
    return {"path": str(poscar_path), "match": bool(dl < 1e-6 and df < 1e-6),
            "max_lattice_diff_A": dl, "max_frac_diff": df}


def formula_units(elements):
    from collections import Counter
    c = Counter(elements)
    g = 0
    for v in c.values():
        g = math.gcd(g, v)
    return g, {k: v // g for k, v in c.items()}


def remove_li(atoms, removed_indices):
    """Return a jarvis Atoms with the given site indices deleted (order kept)."""
    from jarvis.core.atoms import Atoms
    keep = [i for i in range(atoms.num_atoms) if i not in set(removed_indices)]
    return Atoms(lattice_mat=atoms.lattice_mat,
                 elements=[atoms.elements[i] for i in keep],
                 coords=np.array(atoms.frac_coords)[keep], cartesian=False)


def mask_to_removed(mask: int, li_indices):
    return [li_indices[j] for j in range(len(li_indices)) if mask >> j & 1]


# ----------------------------------------------------------------------------
# Worker process: one ALIGNN-FF calculator per process
# ----------------------------------------------------------------------------
_CALC = None
_SETTINGS = None


def _init_worker(settings):
    global _CALC, _SETTINGS
    import torch
    torch.set_num_threads(1)
    from alignn.ff.ff import AlignnAtomwiseCalculator, wt01_path
    _CALC = AlignnAtomwiseCalculator(path=wt01_path(), stress_wt=0.3)   # as dft_prep._alignn_energy
    _SETTINGS = settings


def _sp_energy(jatoms) -> float:
    a = jatoms.ase_converter()
    a.calc = _CALC
    return float(a.get_potential_energy())


def _relax(jatoms, fmax, steps):
    """FIRE on ExpCellFilter (ions + cell).  Returns dict + relaxed jarvis Atoms."""
    from ase.optimize import FIRE
    from ase.constraints import ExpCellFilter
    from jarvis.core.atoms import Atoms
    a = jatoms.ase_converter()
    a.calc = _CALC
    e0 = float(a.get_potential_energy())
    v0 = float(a.get_volume())
    ecf = ExpCellFilter(a)
    opt = FIRE(ecf, logfile=None)
    t = time.time()
    opt.run(fmax=fmax, steps=steps)
    e1 = float(a.get_potential_energy())
    f = np.atleast_2d(a.get_forces())
    fmax_atoms = float(np.linalg.norm(f, axis=1).max()) if len(f) else 0.0
    fmax_filter = float(np.linalg.norm(np.atleast_2d(ecf.get_forces()), axis=1).max())
    relaxed = Atoms(lattice_mat=a.cell.array.tolist(), elements=list(a.get_chemical_symbols()),
                    coords=a.get_scaled_positions().tolist(), cartesian=False)
    return {
        "e_relax": e1, "e_sp_start": e0, "relax_steps": int(opt.get_number_of_steps()),
        "relax_converged": bool(fmax_filter <= fmax), "relax_fmax_atoms": fmax_atoms,
        "relax_fmax_filter": fmax_filter, "relax_volume": float(a.get_volume()),
        "relax_volume_ratio": float(a.get_volume() / v0), "relax_wall_s": time.time() - t,
        "relax_s_per_step": (time.time() - t) / max(1, opt.get_number_of_steps()),
    }, relaxed


def _atoms_dict(j):
    return {"lattice_mat": np.round(np.array(j.lattice_mat), 6).tolist(),
            "elements": list(j.elements),
            "frac_coords": np.round(np.array(j.frac_coords), 6).tolist()}


def _task(args):
    """Dispatch.  args = (kind, payload)."""
    kind, payload = args
    from jarvis.core.atoms import Atoms
    s = _SETTINGS
    base = Atoms.from_dict(s["base_atoms"])
    li = s["li_indices"]
    if kind == "sp":
        mask = payload
        t = time.time()
        e = _sp_energy(remove_li(base, mask_to_removed(mask, li)))
        return kind, {"mask": mask, "e_sp": e, "sp_wall_s": time.time() - t}
    if kind == "relax":
        mask = payload
        info, relaxed = _relax(remove_li(base, mask_to_removed(mask, li)), s["fmax"], s["steps"])
        info["mask"] = mask
        info["relaxed_structure"] = _atoms_dict(relaxed)
        return kind, info
    if kind == "li_ref":
        li_atoms = Atoms.from_dict(payload)
        e_sp = _sp_energy(li_atoms) / li_atoms.num_atoms
        info, relaxed = _relax(li_atoms, s["fmax"], s["steps"])
        return kind, {"structure": "bcc Li, JVASP-913 (JARVIS OptB88-vdW geometry), 2x2x2 supercell",
                      "n_atoms": li_atoms.num_atoms,
                      "e_per_atom_sp": e_sp,
                      "e_per_atom_relaxed": info["e_relax"] / li_atoms.num_atoms,
                      "relax_steps": info["relax_steps"], "relax_converged": info["relax_converged"],
                      "relax_volume_ratio": info["relax_volume_ratio"]}
    if kind == "probe":
        # Sanity of the force field on the pristine cell: analytic forces vs central finite
        # differences, energy jump for small displacements of Li site 0, and the
        # neighbour-shell spacing at the max_neighbors=12 graph cutoff.
        a = base.ase_converter(); a.calc = _CALC
        e0 = float(a.get_potential_energy())
        f = np.atleast_2d(a.get_forces())
        out = {"e_sp": e0, "fmax_eV_A": float(np.linalg.norm(f, axis=1).max()),
               "mean_force_eV_A": float(np.linalg.norm(f, axis=1).mean()), "fd": [], "dE_meV_vs_dx_A": {}}
        probe_atoms = [li[0]] + [k for k in range(base.num_atoms) if base.elements[k] != "Li"][:2]
        for iat in probe_atoms:
            for h in (0.002, 0.01):
                fd = []
                for ax in range(3):
                    es = []
                    for sgn in (1, -1):
                        b = a.copy(); b.calc = _CALC
                        pos = b.get_positions(); pos[iat, ax] += sgn * h; b.set_positions(pos)
                        es.append(float(b.get_potential_energy()))
                    fd.append(-(es[0] - es[1]) / (2 * h))
                out["fd"].append({"atom": iat, "element": base.elements[iat], "h_A": h,
                                  "analytic": f[iat].tolist(), "finite_difference": fd})
        for dx in (0.005, 0.01, 0.02, 0.05, 0.1):
            b = a.copy(); b.calc = _CALC
            pos = b.get_positions(); pos[li[0], 0] += dx; b.set_positions(pos)
            out["dE_meV_vs_dx_A"][str(dx)] = (float(b.get_potential_energy()) - e0) * 1000
        from ase.neighborlist import neighbor_list
        ii, jj, d = neighbor_list("ijd", a, 8.0)
        for iat in probe_atoms:
            ds = np.sort(d[ii == iat])
            out[f"neighbor_distances_atom{iat}"] = {"element": base.elements[iat], "d12_A": float(ds[11]), "d13_A": float(ds[12]),
                                                     "first16_A": np.round(ds[:16], 3).tolist()}
        return kind, out
    if kind == "greedy_contcar":
        # Faithful production loop: rank single points on the previous *relaxed*
        # structure, relax the winner, repeat.  Sequential by construction.
        t0 = time.time()
        cur = base
        removed_sup = []           # supercell indices removed so far (original numbering)
        path = []
        info, cur = _relax(cur, s["fmax"], s["steps"])
        path.append({"step": 0, "n_li": len(li), "removed_li_index": None, "removed_set": [],
                     "mask": 0, "candidate_energies": None, "e_relax": info["e_relax"],
                     "relax_steps": info["relax_steps"], "relax_converged": info["relax_converged"],
                     "relax_volume_ratio": info["relax_volume_ratio"]})
        for step in range(1, len(li) + 1):
            # current structure's Li sites, mapped back to original supercell indices
            remaining = [i for i in li if i not in removed_sup]
            cur_li_pos = [k for k, el in enumerate(cur.elements) if el == "Li"]
            assert len(cur_li_pos) == len(remaining)
            cands = []
            for pos, orig in zip(cur_li_pos, remaining):
                cands.append((orig, pos, _sp_energy(cur.remove_site_by_index(pos))))
            energies = [c[2] for c in cands]
            best = energies.index(min(energies))          # same tie-break as get_next_vacancy
            orig, pos, e_best = cands[best]
            removed_sup.append(orig)
            info, cur = _relax(cur.remove_site_by_index(pos), s["fmax"], s["steps"])
            mask = sum(1 << li.index(i) for i in removed_sup)
            path.append({"step": step, "n_li": len(li) - step, "removed_li_index": orig,
                         "removed_set": sorted(removed_sup), "mask": mask,
                         "candidate_energies": {str(c[0]): c[2] for c in cands},
                         "e_sp_selected_on_prev_relaxed": e_best,
                         "e_relax": info["e_relax"], "relax_steps": info["relax_steps"],
                         "relax_converged": info["relax_converged"],
                         "relax_volume_ratio": info["relax_volume_ratio"]})
        return kind, {"path": path, "wall_s": time.time() - t0}
    raise ValueError(kind)


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------

def greedy_path_from_table(energy_of_mask, n_sites):
    """Greedy on a full energy table: at each step add the single Li whose removal
    gives the lowest energy; ties -> lowest site index (list.index(min) semantics)."""
    mask = 0
    path = [{"step": 0, "n_li": n_sites, "mask": 0, "removed_bit": None,
             "candidate_energies": None, "energy": energy_of_mask[0]}]
    for step in range(1, n_sites + 1):
        cands = [(j, energy_of_mask[mask | (1 << j)]) for j in range(n_sites) if not mask >> j & 1]
        es = [c[1] for c in cands]
        j, e = cands[es.index(min(es))]
        mask |= 1 << j
        path.append({"step": step, "n_li": n_sites - step, "mask": mask, "removed_bit": j,
                     "candidate_energies": {str(c[0]): c[1] for c in cands}, "energy": e})
    return path


def popcount(m):
    return bin(m).count("1")


def per_composition(energy_of_mask, n_sites, path_masks_by_nli, n_fu):
    """Per n_Li statistics; path_masks_by_nli: {name: {n_li: mask}}."""
    out = []
    for n_li in range(n_sites, -1, -1):
        masks = [m for m in energy_of_mask if popcount(m) == n_sites - n_li]
        es = np.array([energy_of_mask[m] for m in masks])
        i_min = int(np.argmin(es))
        e_min = float(es[i_min])
        row = {"n_li": n_li, "n_configs": len(masks), "min_mask": masks[i_min], "e_min": e_min,
               "n_within_25meV": int(np.sum(es - e_min <= GAP_WINDOW_MEV / 1000.0)),
               "n_distinct_levels": int(len(np.unique(np.round(es / DEGENERACY_TOL_EV)))),
               "n_degenerate_with_min": int(np.sum(np.abs(es - e_min) <= DEGENERACY_TOL_EV)),
               "spread_meV": float((es.max() - e_min) * 1000)}
        for name, bym in path_masks_by_nli.items():
            m = bym.get(n_li)
            if m is None or m not in energy_of_mask:
                row[name] = None
                continue
            e = energy_of_mask[m]
            gap = (e - e_min) * 1000.0
            rank = int(np.sum(es < e - 1e-12)) + 1
            row[name] = {"mask": m, "energy": e, "gap_meV_cell": gap, "gap_meV_fu": gap / n_fu,
                         "found_min": bool(gap <= 1.0),        # 1 meV = float32 noise floor
                         "exact_same_config": bool(m == masks[i_min]), "rank": rank}
        out.append(row)
    return out


def lower_hull_voltages(points, e_li):
    """points: list of (n_li, E) with one entry per n_li, n_li descending.
    Returns step voltages and convex-hull plateau voltages (dft_prep convention:
    V = (E_lo - E_hi)/dn + e_li)."""
    pts = sorted(points, key=lambda p: p[0])            # n ascending
    n_tot = pts[-1][0]
    e_empty, e_full = pts[0][1], pts[-1][1]
    fe = [(n, E - (n / n_tot) * e_full - (1 - n / n_tot) * e_empty, E) for n, E in pts]
    hull = []
    for p in fe:                                        # monotone chain, lower hull
        while len(hull) >= 2:
            (x0, y0, _), (x1, y1, _) = hull[-2], hull[-1]
            x2, y2, _ = p
            if (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0) <= 1e-12:
                hull.pop()
            else:
                break
        hull.append(p)
    plateaus = []
    for (n_lo, _, e_lo), (n_hi, _, e_hi) in zip(hull[:-1], hull[1:]):
        plateaus.append({"n_li_from": n_hi, "n_li_to": n_lo, "x_from": n_hi / n_tot, "x_to": n_lo / n_tot,
                         "voltage": (e_lo - e_hi) / (n_hi - n_lo) + e_li})
    plateaus.sort(key=lambda p: -p["n_li_from"])
    steps = []
    for (n_lo, E_lo), (n_hi, E_hi) in zip(pts[:-1], pts[1:]):
        steps.append({"n_li_from": n_hi, "n_li_to": n_lo, "voltage": (E_lo - E_hi) / (n_hi - n_lo) + e_li})
    steps.sort(key=lambda p: -p["n_li_from"])
    avg = (e_empty - e_full) / n_tot + e_li
    return {"e_li_ref": e_li, "average_voltage": avg, "step_voltages": steps,
            "hull_plateaus": plateaus, "formation_energies_meV_cell": {str(n): y * 1000 for n, y, _ in fe},
            "hull_vertices_n_li": [n for n, _, _ in hull]}


def subset_chain_checks(path_masks_by_nli, min_masks_by_nli, n_sites):
    """Is greedy config at step k (n_li) a subset of the global-min config at step k+1
    (n_li - 1)?  Removed-set inclusion."""
    rows = []
    for n_li in range(n_sites, 0, -1):
        g = path_masks_by_nli.get(n_li)
        m_next = min_masks_by_nli.get(n_li - 1)
        m_here = min_masks_by_nli.get(n_li)
        if g is None or m_next is None:
            continue
        rows.append({"n_li": n_li, "greedy_mask": g, "globalmin_next_mask": m_next,
                     "greedy_subset_of_next_min": bool(g & m_next == g),
                     "greedy_equals_min_here": bool(g == m_here),
                     "min_here_subset_of_next_min": bool(m_here & m_next == m_here)})
    return rows


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def run(jid, args):
    import multiprocessing as mp
    from jarvis.core.atoms import Atoms
    from jarvis.analysis.thermodynamics.energetics import unary_energy

    out_json = HERE / f"results_{jid}.json"
    log = lambda *a: print(time.strftime("%H:%M:%S"), *a, flush=True)

    sup, dim, row = build_supercell(jid, args.cache)
    li_indices = [i for i, el in enumerate(sup.elements) if el == "Li"]
    n_sites = len(li_indices)
    n_fu, fu = formula_units(sup.elements)
    ref_p = Path(args.ref_poscar) if args.ref_poscar else find_reference_poscar(jid)
    ref_cmp = compare_to_poscar(sup, ref_p) if ref_p else None
    log(f"{jid} {row['formula']} spg {row['spg_number']} dim {dim} atoms {sup.num_atoms} Li {n_sites} "
        f"f.u./cell {n_fu} ref_match {ref_cmp}")
    if n_sites > 12:
        raise SystemExit(f"{n_sites} Li sites -> {2**n_sites} configs; refusing (>12).")

    li_prim, _ = load_from_cache("JVASP-913", args.cache)
    li_atoms = li_prim.make_supercell([2, 2, 2])      # 8-atom bcc cell (1-atom cell trips the calculator's force shape)
    settings = {"base_atoms": sup.to_dict(), "li_indices": li_indices,
                "fmax": args.fmax, "steps": args.steps}

    result = {}
    if out_json.exists() and args.resume:
        result = json.load(open(out_json))
        log(f"resuming from {out_json}: {len(result.get('configs', {}))} configs present")
    configs = {int(k): v for k, v in result.get("configs", {}).items()}
    for m in range(2 ** n_sites):
        configs.setdefault(m, {"mask": m, "removed_li_indices": mask_to_removed(m, li_indices),
                               "n_li": n_sites - popcount(m)})

    def save():
        result.update({
            "jid": jid, "formula": row["formula"], "spg_number": row["spg_number"], "dim": dim,
            "n_atoms": sup.num_atoms, "n_li": n_sites, "n_configs": 2 ** n_sites,
            "formula_units_per_cell": n_fu, "formula_unit": fu, "li_indices": li_indices,
            "li_frac_coords": np.round(np.array(sup.frac_coords)[li_indices], 6).tolist(),
            "reference_poscar": ref_cmp,
            "settings": {"alignn_weights": "alignnff_wt01 (wt01_path)", "stress_wt": 0.3,
                         "relax": "ASE FIRE on ExpCellFilter (ions + cell)",
                         "fmax_eV_per_A": args.fmax, "step_cap": args.steps,
                         "relax_mode_requested": args.relax, "workers": args.workers,
                         "gap_window_meV": GAP_WINDOW_MEV, "degeneracy_tol_eV": DEGENERACY_TOL_EV},
            "e_li_reference": result.get("e_li_reference", {}),
            "force_probe": result.get("force_probe"),
            "configs": {str(m): configs[m] for m in sorted(configs)},
        })
        tmp = out_json.with_suffix(".json.tmp")
        json.dump(result, open(tmp, "w"), indent=1)
        os.replace(tmp, out_json)

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(args.workers, initializer=_init_worker, initargs=(settings,))
    timing = result.setdefault("timing", {})

    # ---- Li reference + level (a): all single points ---------------------------
    todo = [("sp", m) for m in sorted(configs) if "e_sp" not in configs[m]]
    if "e_li_reference" not in result or "alignn_ff_bcc_li" not in result["e_li_reference"]:
        todo = [("li_ref", li_atoms.to_dict())] + todo
    if "force_probe" not in result:
        todo = [("probe", None)] + todo
    log(f"level (a): {len(todo)} tasks on {args.workers} workers")
    t0 = time.time()
    for k, (kind, payload) in enumerate(pool.imap_unordered(_task, todo), 1):
        if kind == "li_ref":
            result["e_li_reference"] = {
                "alignn_ff_bcc_li": payload,
                "used_level_a_sp": payload["e_per_atom_sp"],
                "used_level_b_relaxed": payload["e_per_atom_relaxed"],
                "screen_py_unary_energy_Li": float(unary_energy("Li")),
                "dft_prep_E_LI_METAL": E_LI_DFT,
                "note": ("Voltages in this file use the ALIGNN-FF energy of bcc Li computed with the same "
                         "calculator: single-point on the JVASP-913 geometry for level (a), FIRE-relaxed "
                         "for level (b).  screen.py uses jarvis unary_energy('Li') (JARVIS OptB88-vdW DFT); "
                         "dft_prep.py uses in-house Li_sv DFT values.  The choice shifts every voltage by a "
                         "constant and does not affect any energy gap or ranking result."),
            }
            log(f"Li ref: sp {payload['e_per_atom_sp']:.4f}  relaxed {payload['e_per_atom_relaxed']:.4f} eV/atom")
        elif kind == "probe":
            result["force_probe"] = payload
            log(f"probe: pristine fmax {payload['fmax_eV_A']:.1f} eV/A, dE(0.01 A on Li0) {payload['dE_meV_vs_dx_A']['0.01']:.1f} meV")
        else:
            configs[payload["mask"]].update(payload)
        if k % 32 == 0 or k == len(todo):
            log(f"  {k}/{len(todo)} done, {time.time()-t0:.0f}s")
            save()
    timing["level_a_wall_s"] = timing.get("level_a_wall_s", 0.0) + time.time() - t0
    sp = {m: configs[m]["e_sp"] for m in configs}
    timing["sp_mean_s"] = float(np.mean([configs[m]["sp_wall_s"] for m in configs if "sp_wall_s" in configs[m]] or [0]))

    # ---- greedy on single points (exact get_next_vacancy on the rigid cell) -----
    gsp = greedy_path_from_table(sp, n_sites)
    for p in gsp:
        p["removed_set"] = mask_to_removed(p["mask"], li_indices)
        p["removed_li_index"] = None if p["removed_bit"] is None else li_indices[p["removed_bit"]]
    result["greedy_sp_path"] = gsp
    save()

    # ---- level (b): relaxations ----------------------------------------------------
    relax_mode = args.relax
    minimal_masks = set()
    comp_sp = per_composition(sp, n_sites, {"greedy_sp": {p["n_li"]: p["mask"] for p in gsp}}, n_fu)
    for r in comp_sp:
        minimal_masks.add(r["min_mask"])
        minimal_masks.add(r["greedy_sp"]["mask"])
    if relax_mode != "none":
        # greedy_contcar runs as one sequential task alongside the enumeration
        gc_needed = "greedy_contcar_path" not in result
        pending = [m for m in sorted(configs) if "e_relax" not in configs[m]]
        if relax_mode == "minimal":
            pending = [m for m in pending if m in minimal_masks]
        # Put the minimal set first so a budget abort still leaves the useful subset done,
        # then the first-10 timing probe is on ordinary configs.
        pending.sort(key=lambda m: (m not in minimal_masks, m))
        tasks = ([("greedy_contcar", None)] if gc_needed else []) + [("relax", m) for m in pending]
        log(f"level (b) mode={relax_mode}: {len(pending)} relaxations (+greedy_contcar={gc_needed})")
        t0 = time.time()
        n_done = 0
        probe_reported = False
        decided_minimal = relax_mode == "minimal"
        it = pool.imap_unordered(_task, tasks)
        walls = []
        for kind, payload in it:
            if kind == "greedy_contcar":
                result["greedy_contcar_path"] = payload["path"]
                timing["greedy_contcar_wall_s"] = payload["wall_s"]
                log(f"greedy_contcar done in {payload['wall_s']:.0f}s")
                save()
                continue
            m = payload.pop("mask")
            configs[m].update(payload)
            n_done += 1
            walls.append(payload["relax_wall_s"])
            if n_done % 8 == 0 or n_done == len(pending):
                log(f"  relax {n_done}/{len(pending)} elapsed {time.time()-t0:.0f}s "
                    f"mean/config {np.mean(walls):.0f}s steps~{payload['relax_steps']}")
                save()
            if n_done == 10 and not probe_reported:
                probe_reported = True
                elapsed = time.time() - t0
                proj = elapsed / 10 * 2 ** n_sites
                timing["relax_probe_first10_wall_s"] = elapsed
                timing["relax_projected_all256_wall_s"] = proj
                log(f"PROBE: first 10 relaxations took {elapsed:.0f}s wall on {args.workers} workers "
                    f"-> projected {proj/3600:.2f} h for all {2**n_sites} (budget {args.budget_hours} h)")
                if relax_mode == "auto" and proj > args.budget_hours * 3600:
                    decided_minimal = True
                    log("  over budget -> relaxing only the per-composition SP-minimum and greedy configs")
                    # Drain: cannot cancel imap tasks individually; terminate and resubmit the minimal rest.
                    pool.terminate(); pool.join()
                    pool = ctx.Pool(args.workers, initializer=_init_worker, initargs=(settings,))
                    rest = [m for m in pending if m in minimal_masks and "e_relax" not in configs[m]]
                    tasks2 = ([("greedy_contcar", None)] if "greedy_contcar_path" not in result else []) + \
                             [("relax", m) for m in rest]
                    pending = [m for m in pending if "e_relax" in configs[m]] + rest
                    it = pool.imap_unordered(_task, tasks2)
                    for kind2, payload2 in it:
                        if kind2 == "greedy_contcar":
                            result["greedy_contcar_path"] = payload2["path"]
                            timing["greedy_contcar_wall_s"] = payload2["wall_s"]
                            save(); continue
                        m2 = payload2.pop("mask"); configs[m2].update(payload2); n_done += 1
                        log(f"  relax(minimal) {n_done}/{len(pending)}")
                        save()
                    break
        timing["level_b_wall_s"] = timing.get("level_b_wall_s", 0.0) + time.time() - t0
        result["settings"]["relax_mode_effective"] = "minimal" if decided_minimal else "all"
    pool.terminate(); pool.join()

    # ---- analysis ---------------------------------------------------------------------
    e_li_a = result["e_li_reference"]["used_level_a_sp"]
    e_li_b = result["e_li_reference"]["used_level_b_relaxed"]
    paths_a = {"greedy_sp": {p["n_li"]: p["mask"] for p in gsp}}
    if "greedy_contcar_path" in result:
        paths_a["greedy_contcar"] = {p["n_li"]: p["mask"] for p in result["greedy_contcar_path"]}
    comp_a = per_composition(sp, n_sites, paths_a, n_fu)
    min_a = {r["n_li"]: r["min_mask"] for r in comp_a}
    analysis = {"level_a_single_point": {
        "per_composition": comp_a,
        "greedy_sp_found_min_all": all(r["greedy_sp"]["found_min"] for r in comp_a),
        "subset_chain": subset_chain_checks(paths_a["greedy_sp"], min_a, n_sites),
        "voltage": {
            "greedy_sp_path": lower_hull_voltages([(p["n_li"], sp[p["mask"]]) for p in gsp], e_li_a),
            "global_min_path": lower_hull_voltages([(r["n_li"], r["e_min"]) for r in comp_a], e_li_a),
            "e_li_ref_alt_screen_py": E_LI_UNARY_SCREEN_PY,
        }}}
    rel = {m: configs[m]["e_relax"] for m in configs if "e_relax" in configs[m]}
    if rel:
        full = len(rel) == 2 ** n_sites
        paths_b = dict(paths_a)
        comp_b = per_composition(rel, n_sites, paths_b, n_fu)
        min_b = {r["n_li"]: r["min_mask"] for r in comp_b}
        lvl = {"complete_enumeration": full, "n_relaxed": len(rel), "per_composition": comp_b,
               "n_not_converged": int(sum(1 for m in rel if not configs[m]["relax_converged"])),
               "max_volume_ratio": float(max(configs[m]["relax_volume_ratio"] for m in rel)),
               "min_volume_ratio": float(min(configs[m]["relax_volume_ratio"] for m in rel)),
               "subset_chain_greedy_sp": subset_chain_checks(paths_a["greedy_sp"], min_b, n_sites),
               "voltage": {"global_min_relaxed_path": lower_hull_voltages([(r["n_li"], r["e_min"]) for r in comp_b], e_li_b)}}
        if all(paths_a["greedy_sp"][n] in rel for n in paths_a["greedy_sp"]):
            lvl["voltage"]["greedy_sp_path_relaxed"] = lower_hull_voltages(
                [(n, rel[m]) for n, m in paths_a["greedy_sp"].items()], e_li_b)
        if "greedy_contcar_path" in result:
            gc = result["greedy_contcar_path"]
            lvl["subset_chain_greedy_contcar"] = subset_chain_checks(paths_a["greedy_contcar"], min_b, n_sites)
            lvl["voltage"]["greedy_contcar_path"] = lower_hull_voltages([(p["n_li"], p["e_relax"]) for p in gc], e_li_b)
            # greedy_contcar's own relaxed energies vs enumeration minima (its structures come from a
            # chained relaxation, so they can differ slightly from the fresh relax of the same mask)
            for r in comp_b:
                p = next(p for p in gc if p["n_li"] == r["n_li"])
                gap = (p["e_relax"] - r["e_min"]) * 1000
                r["greedy_contcar_chained"] = {"mask": p["mask"], "energy": p["e_relax"], "gap_meV_cell": gap,
                                               "gap_meV_fu": gap / n_fu, "found_min": bool(gap <= 1.0),
                                               "exact_same_config": bool(p["mask"] == r["min_mask"])}
            lvl["greedy_contcar_found_min_all"] = all(r["greedy_contcar_chained"]["found_min"] for r in comp_b)
        if full:
            gr = greedy_path_from_table(rel, n_sites)
            lvl["greedy_on_relaxed_energies_path"] = gr
            paths_b["greedy_relaxed"] = {p["n_li"]: p["mask"] for p in gr}
            lvl["per_composition"] = per_composition(rel, n_sites, paths_b, n_fu)
            lvl["greedy_sp_found_min_all"] = all(r["greedy_sp"]["found_min"] for r in lvl["per_composition"])
        analysis["level_b_relaxed"] = lvl
    result["analysis"] = analysis
    save()
    log(f"wrote {out_json}")
    return result


# ----------------------------------------------------------------------------
# Plot + summary
# ----------------------------------------------------------------------------
C_ALL, C_GREEDY, C_MIN, C_CONTCAR, C_HULL, C_TEXT = "#b5b4ad", "#eb6834", "#2a78d6", "#4a3aa7", "#0b0b0b", "#52514e"


def _fe_curve(energy_by_mask_or_list, e_full, e_empty, n_tot):
    def fe(n, E):
        return (E - (n / n_tot) * e_full - (1 - n / n_tot) * e_empty) * 1000
    return fe


def plot_system(res, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 15,
                         "legend.fontsize": 13, "xtick.labelsize": 14, "ytick.labelsize": 14})
    n = res["n_li"]
    cfg = {int(k): v for k, v in res["configs"].items()}
    a_sp = res["analysis"]["level_a_single_point"]
    have_b = "level_b_relaxed" in res["analysis"]
    ncol = 3 if have_b else 2
    fig, axes = plt.subplots(1, ncol, figsize=(7.6 * ncol, 6.4))
    gsp = res["greedy_sp_path"]

    def mixing_panel(ax, key, comp, title, vkey, contcar=False):
        have = {m: v[key] for m, v in cfg.items() if key in v}
        e_min = {r["n_li"]: r["e_min"] for r in comp}
        e_full, e_empty = e_min[n], e_min[0]
        fe = lambda x, E: (E - (x / n) * e_full - (1 - x / n) * e_empty) * 1000
        xs = [n - bin(m).count("1") for m in have]
        ax.scatter(xs, [fe(x, have[m]) for x, m in zip(xs, have)], s=22, c=C_ALL, zorder=1,
                   label=f"all {len(have)} configurations" if len(have) == 2 ** n else f"{len(have)} relaxed configurations")
        gx = list(range(n, -1, -1))
        ax.plot(gx, [fe(x, e_min[x]) for x in gx], "-", color=C_MIN, lw=2, marker="o", ms=8, zorder=3,
                label="global minimum per n$_{Li}$")
        py = [(p["n_li"], fe(p["n_li"], have[p["mask"]])) for p in gsp if p["mask"] in have]
        ax.plot([q[0] for q in py], [q[1] for q in py], "--", color=C_GREEDY, lw=2, marker="s", ms=8, mfc="none", mew=2,
                zorder=4, label="greedy (single-point ranking)")
        if contcar and "greedy_contcar_path" in res:
            gc = res["greedy_contcar_path"]
            ax.plot([p["n_li"] for p in gc], [fe(p["n_li"], p["e_relax"]) for p in gc], ":", color=C_CONTCAR,
                    lw=2.2, marker="^", ms=8, mfc="none", mew=2, zorder=5, label="greedy (relax-then-rank loop)")
        hv = res["analysis"][vkey[0]]["voltage"][vkey[1]]["hull_vertices_n_li"]
        ax.plot(hv, [fe(x, e_min[x]) for x in hv], "-", color=C_HULL, lw=1, alpha=0.6, zorder=2, label="convex hull")
        ax.set_xlabel("n$_{Li}$ remaining in cell"); ax.set_ylabel("mixing energy (meV / cell)")
        ax.set_title(title); ax.set_xticks(range(0, n + 1)); ax.invert_xaxis(); ax.grid(alpha=0.25)
        ax.legend(loc="upper center", frameon=False)

    mixing_panel(axes[0], "e_sp", a_sp["per_composition"], "(a) single-point, unrelaxed lattice",
                 ("level_a_single_point", "global_min_path"))
    # gap panel: E - E_min(n_Li), symlog so 25 meV and 1 eV are both visible
    ax = axes[1]
    e_min = {r["n_li"]: r["e_min"] for r in a_sp["per_composition"]}
    xs = [n - bin(m).count("1") for m in cfg]
    ys = [(cfg[m]["e_sp"] - e_min[n - bin(m).count("1")]) * 1000 for m in cfg]
    ax.scatter(xs, ys, s=22, c=C_ALL, zorder=1, label="all configurations")
    gy = [(cfg[p["mask"]]["e_sp"] - e_min[p["n_li"]]) * 1000 for p in gsp]
    ax.plot([p["n_li"] for p in gsp], gy, "--", color=C_GREEDY, lw=2, marker="s", ms=9, mfc="none", mew=2, zorder=4,
            label="greedy configuration")
    for p, y in zip(gsp, gy):
        if y > 1.0:
            ax.annotate(f"{y:.0f}", (p["n_li"], y), textcoords="offset points", xytext=(8, 4), fontsize=13, color=C_TEXT)
    ax.axhline(GAP_WINDOW_MEV, color=C_MIN, lw=1.2, ls=":", label=f"{GAP_WINDOW_MEV:.0f} meV window")
    ax.set_yscale("symlog", linthresh=GAP_WINDOW_MEV, linscale=1.0)
    ax.set_ylim(-3, max(ys) * 1.6)
    ax.set_xlabel("n$_{Li}$ remaining in cell"); ax.set_ylabel("E $-$ E$_{min}$(n$_{Li}$)  (meV / cell, symlog)")
    ax.set_title("(a') greedy gap to the exhaustive minimum"); ax.set_xticks(range(0, n + 1)); ax.invert_xaxis()
    ax.grid(alpha=0.25, which="both"); ax.legend(loc="upper center", frameon=False)
    if have_b:
        mixing_panel(axes[2], "e_relax", res["analysis"]["level_b_relaxed"]["per_composition"],
                     "(b) after FIRE steps (diagnostic, see caveat)", ("level_b_relaxed", "global_min_relaxed_path"), contcar=True)
    fig.suptitle(f"{res['jid']}  {res['formula']}  ({res['n_atoms']} atoms, {n} Li sites, {2**n} configurations)", fontsize=15)
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def _fmt_v(v):
    return "; ".join(f"{p['voltage']:.3f} V (x {p['x_from']:.3f}→{p['x_to']:.3f})" for p in v["hull_plateaus"])


def summarize(jids):
    lines = ["# Exhaustive Li-vacancy enumeration vs. greedy ALIGNN-FF ranking (R1.M6)", ""]
    lines += ["Generated by `enumerate.py`.  For each cell every subset of the Li sites (2^8 = 256 "
              "configurations, no symmetry reduction) was evaluated with ALIGNN-FF (`alignnff_wt01`, "
              "same calculator as `dft_prep._alignn_energy`).  `greedy_sp` is the `get_next_vacancy` "
              "heuristic applied to the rigid step_00 cell; `greedy_contcar` is the production loop "
              "(rank single points on the previous step's relaxed structure, relax the winner, repeat) "
              "with ALIGNN-FF standing in for VASP.  Gaps are energy of the greedy configuration minus "
              "the exhaustive minimum at the same n_Li; `found` means gap <= 1 meV/cell.  Voltages use the "
              "ALIGNN-FF bcc-Li reference computed with the same calculator (values listed per system); "
              "screen.py instead uses jarvis `unary_energy('Li') = -0.925` eV/atom (JARVIS OptB88-vdW DFT) "
              "and dft_prep.py uses in-house Li_sv DFT values; a different reference shifts all voltages by "
              "a constant and changes no gap or ranking.", ""]
    for jid in jids:
        p = HERE / f"results_{jid}.json"
        if not p.exists():
            lines += [f"## {jid}: results_{jid}.json missing", ""]
            continue
        r = json.load(open(p))
        if "analysis" not in r:
            lines += [f"## {jid}: results_{jid}.json still in progress (no analysis block yet)", ""]
            continue
        a = r["analysis"]
        n, nfu = r["n_li"], r["formula_units_per_cell"]
        lines += [f"## {jid}  {r['formula']}", ""]
        lines += [f"- cell: {r['dim'][0]}x{r['dim'][1]}x{r['dim'][2]} supercell, {r['n_atoms']} atoms, {n} Li sites, "
                  f"{nfu} formula units per cell; reference POSCAR match: "
                  f"{r['reference_poscar']['match'] if r['reference_poscar'] else 'n/a'} ({r['reference_poscar']['path'] if r['reference_poscar'] else ''})"]
        eli = r["e_li_reference"]
        lines += [f"- Li reference (ALIGNN-FF, bcc Li JVASP-913): single-point {eli['used_level_a_sp']:.4f} eV/atom "
                  f"(level a), relaxed {eli['used_level_b_relaxed']:.4f} eV/atom (level b); screen.py unary_energy = "
                  f"{eli['screen_py_unary_energy_Li']:.3f}; dft_prep DFT refs {eli['dft_prep_E_LI_METAL']}"]
        t = r.get("timing", {})
        s = r["settings"]
        lines += [f"- relax settings: {s['relax']}, fmax {s['fmax_eV_per_A']} eV/A, step cap {s['step_cap']}, "
                  f"mode requested `{s['relax_mode_requested']}`, effective `{s.get('relax_mode_effective', 'none')}`; "
                  f"timing: single-point {t.get('sp_mean_s', float('nan')):.2f} s/config, level (a) wall "
                  f"{t.get('level_a_wall_s', 0)/60:.1f} min, first-10 relax probe {t.get('relax_probe_first10_wall_s', float('nan')):.0f} s "
                  f"-> projected {t.get('relax_projected_all256_wall_s', float('nan'))/3600:.2f} h for 256, "
                  f"level (b) wall {t.get('level_b_wall_s', 0)/3600:.2f} h on {s['workers']} workers", ""]
        fp = r.get("force_probe")
        if fp:
            nb = [v for k, v in fp.items() if k.startswith("neighbor_distances_atom")]
            lines += [f"- force-field sanity probe (pristine cell): analytic fmax {fp['fmax_eV_A']:.1f} eV/A, mean |F| {fp['mean_force_eV_A']:.1f} eV/A "
                      f"on a DFT-relaxed structure; energy change for displacing Li site 0 by 0.005/0.01/0.02/0.05/0.1 A: "
                      + "/".join(f"{fp['dE_meV_vs_dx_A'][k]:+.0f}" for k in ("0.005", "0.01", "0.02", "0.05", "0.1")) + " meV/cell; "
                      f"12th/13th-neighbour distances: " + ", ".join(f"{x['element']} {x['d12_A']:.3f}/{x['d13_A']:.3f} A" for x in nb)]
            lines += ["- analytic vs finite-difference forces (eV/A): " + "; ".join(
                f"{x['element']}{x['atom']} h={x['h_A']}: analytic ({', '.join(f'{v:.2f}' for v in x['analytic'])}) FD ({', '.join(f'{v:.2f}' for v in x['finite_difference'])})"
                for x in fp["fd"]), ""]
        # level a
        lines += ["### Level (a): single-point energies, unrelaxed lattice", ""]
        ca = a["level_a_single_point"]["per_composition"]
        c7 = next(c for c in ca if c["n_li"] == n - 1)
        lines += [f"Noise floor: the {c7['n_configs']} single-vacancy configurations at n_Li = {n-1} are all symmetry-equivalent in this cell, "
                  f"yet ALIGNN-FF spreads them over {c7['spread_meV']:.1f} meV/cell ({c7['n_distinct_levels']} distinct energies); "
                  f"gaps below this value are within the model's own symmetry-breaking noise.  'rank' counts configurations with strictly lower "
                  f"energy plus one, so a rank of 2 with a 0.0 meV gap is a degenerate tie.", ""]
        lines += ["| n_Li | configs | E_min (eV) | greedy E (eV) | gap meV/cell | gap meV/f.u. | greedy found min | rank of greedy | within 25 meV | distinct levels | degenerate w/ min | spread max-min meV |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for c in ca:
            g = c["greedy_sp"]
            lines.append(f"| {c['n_li']} | {c['n_configs']} | {c['e_min']:.4f} | {g['energy']:.4f} | {g['gap_meV_cell']:.1f} | "
                         f"{g['gap_meV_fu']:.1f} | {'yes' if g['found_min'] else 'NO'} | {g['rank']}/{c['n_configs']} | "
                         f"{c['n_within_25meV']} | {c['n_distinct_levels']} | {c['n_degenerate_with_min']} | {c['spread_meV']:.1f} |")
        sc = a["level_a_single_point"]["subset_chain"]
        lines += ["", "Path nesting (removed-set inclusion): greedy config at n_Li a subset of the global-min config at n_Li-1?  " +
                  ", ".join(f"{x['n_li']}→{x['n_li']-1}: {'yes' if x['greedy_subset_of_next_min'] else 'no'}" for x in sc)]
        lines += ["Global-min configs nested along a single path (min at n_Li a subset of min at n_Li-1)?  " +
                  ", ".join(f"{x['n_li']}→{x['n_li']-1}: {'yes' if x['min_here_subset_of_next_min'] else 'no'}" for x in sc), ""]
        va = a["level_a_single_point"]["voltage"]
        lines += ["| path | average voltage (V) | hull plateaus |", "|---|---|---|",
                  f"| greedy (single-point ranking) | {va['greedy_sp_path']['average_voltage']:.3f} | {_fmt_v(va['greedy_sp_path'])} |",
                  f"| global minimum per composition | {va['global_min_path']['average_voltage']:.3f} | {_fmt_v(va['global_min_path'])} |", ""]
        lines += ["Step voltages (V) n_Li 8→7 ... 1→0, greedy: " + ", ".join(f"{x['voltage']:.3f}" for x in va['greedy_sp_path']['step_voltages'])]
        lines += ["Step voltages (V) n_Li 8→7 ... 1→0, global-min: " + ", ".join(f"{x['voltage']:.3f}" for x in va['global_min_path']['step_voltages']), ""]
        # level b
        if "level_b_relaxed" in a:
            b = a["level_b_relaxed"]
            lines += ["**Caveat.** The force-field sanity probe above shows that on these cells the calculator's analytic forces are "
                      "10-20x larger than, and inconsistent in sign with, finite-difference forces, and that the energy is discontinuous "
                      "at the 0.01 A scale (degenerate neighbour shells at the max_neighbors=12 graph cutoff).  FIRE therefore never reaches "
                      "fmax = 0.05 eV/A (the residual is ~10 eV/A) and the step cap always binds, so a 'relaxed' energy here is the energy "
                      "after N FIRE steps on a rough surface, not a converged local minimum.  These numbers are reported as a diagnostic only; "
                      "the production workflow never relaxes with ALIGNN-FF (VASP does the relaxation and ALIGNN-FF supplies only the "
                      "single-point ranking of level (a)).", ""]
            lines += [f"### Level (b): ALIGNN-FF relaxed (ions + cell), {b['n_relaxed']}/{2**n} configurations relaxed"
                      f" ({'complete' if b['complete_enumeration'] else 'partial: per-composition SP-minimum + greedy configs only'})", ""]
            lines += [f"- not converged within step cap: {b['n_not_converged']}; relaxed/initial volume ratio range "
                      f"{b['min_volume_ratio']:.3f}–{b['max_volume_ratio']:.3f}", ""]
            hdr = "| n_Li | relaxed configs | E_min relaxed (eV) | greedy_sp config: E (eV) / gap meV/cell / gap meV/f.u. / found / rank | " \
                  "greedy_contcar (chained): E (eV) / gap meV/cell / gap meV/f.u. / found | within 25 meV | distinct levels |"
            lines += [hdr, "|" + "---|" * 7]
            for c in b["per_composition"]:
                g = c.get("greedy_sp"); gc = c.get("greedy_contcar_chained")
                gs = f"{g['energy']:.4f} / {g['gap_meV_cell']:.1f} / {g['gap_meV_fu']:.1f} / {'yes' if g['found_min'] else 'NO'} / {g['rank']}/{c['n_configs']}" if g else "n/a"
                gcs = f"{gc['energy']:.4f} / {gc['gap_meV_cell']:.1f} / {gc['gap_meV_fu']:.1f} / {'yes' if gc['found_min'] else 'NO'}" if gc else "n/a"
                lines.append(f"| {c['n_li']} | {c['n_configs']} | {c['e_min']:.4f} | {gs} | {gcs} | {c['n_within_25meV']} | {c['n_distinct_levels']} |")
            if "greedy_on_relaxed_energies_path" in b:
                lines += ["", "Greedy re-run on the relaxed-energy table (hypothetical: rank by relaxed energy instead of single point): " +
                          ", ".join(f"n_Li {c['n_li']}: gap {c['greedy_relaxed']['gap_meV_cell']:.1f} meV" for c in b["per_composition"] if c.get("greedy_relaxed"))]
            sc = b["subset_chain_greedy_sp"]
            lines += ["", "Nesting vs relaxed global minima, greedy_sp: " + ", ".join(f"{x['n_li']}→{x['n_li']-1}: {'yes' if x['greedy_subset_of_next_min'] else 'no'}" for x in sc)]
            if "subset_chain_greedy_contcar" in b:
                sc = b["subset_chain_greedy_contcar"]
                lines += ["Nesting vs relaxed global minima, greedy_contcar: " + ", ".join(f"{x['n_li']}→{x['n_li']-1}: {'yes' if x['greedy_subset_of_next_min'] else 'no'}" for x in sc)]
            lines += ["Relaxed global-min configs nested along one path?  " + ", ".join(f"{x['n_li']}→{x['n_li']-1}: {'yes' if x['min_here_subset_of_next_min'] else 'no'}" for x in sc), ""]
            vb = b["voltage"]
            lines += ["| path (relaxed energies) | average voltage (V) | hull plateaus |", "|---|---|---|"]
            for key, lab in [("greedy_sp_path_relaxed", "greedy_sp configs, relaxed"), ("greedy_contcar_path", "greedy relax-then-rank loop"),
                             ("global_min_relaxed_path", "global minimum per composition (relaxed)")]:
                if key in vb:
                    lines.append(f"| {lab} | {vb[key]['average_voltage']:.3f} | {_fmt_v(vb[key])} |")
            lines += [""]
            for key, lab in [("greedy_contcar_path", "greedy loop"), ("global_min_relaxed_path", "global-min")]:
                if key in vb:
                    lines += [f"Step voltages (V) 8→7 ... 1→0, {lab}: " + ", ".join(f"{x['voltage']:.3f}" for x in vb[key]['step_voltages'])]
            lines += [""]
        lines += [f"![enumeration](enumeration_{jid}.png)", ""]
    (HERE / "enumeration_summary.md").write_text("\n".join(lines))
    print("wrote", HERE / "enumeration_summary.md")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jid", nargs="?", help="JARVIS ID, e.g. JVASP-2017")
    ap.add_argument("--relax", default="auto", choices=["auto", "all", "minimal", "none"])
    ap.add_argument("--budget-hours", type=float, default=3.0, help="per-system wall budget for --relax auto")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--fmax", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    ap.add_argument("--ref-poscar", default=None)
    ap.add_argument("--resume", action="store_true", help="reuse energies already in results_JID.json")
    ap.add_argument("--summary", action="store_true", help="only rebuild summary + plots from existing JSON")
    ap.add_argument("--jids", nargs="*", default=["JVASP-2017", "JVASP-96563"])
    args = ap.parse_args()
    if args.summary:
        for jid in args.jids:
            p = HERE / f"results_{jid}.json"
            if p.exists():
                plot_system(json.load(open(p)), HERE / f"enumeration_{jid}.png")
        summarize(args.jids)
        return
    if not args.jid:
        ap.error("jid required unless --summary")
    res = run(args.jid, args)
    plot_system(res, HERE / f"enumeration_{args.jid}.png")
    summarize(args.jids)


if __name__ == "__main__":
    main()
