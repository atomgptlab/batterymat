"""Voltages, hull plateaus and plots for every revision DFT campaign.

For each campaign directory under ../dft_inputs/ that has at least two recorded
energies this script

  - regenerates the discharge-curve and voltage-curve plots
    (dft_prep.compute_voltage_curve, saved in this folder as
    discharge_curve_<id>.png / voltage_curve_<id>.png),
  - computes the lower-convex-hull plateaus over the recorded steps,
  - computes the average voltage over the full range and over the range that
    matches the experimental comparison (the average over any range telescopes
    to the two-endpoint value (E_lo - E_hi)/(z dn) + mu/z),
  - flags steps whose results/OUTCAR lacks "reached required accuracy",

and writes revision_voltages.json and revision_voltages.md. Metal references
come from energies.json, falling back to dft_prep._E_ION_METAL.

Usage (from screening_cathode/analysis/):  python revision_voltages.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import dft_prep as dp  # noqa: E402

# dir name, label, experimental comparison range as (n_hi, n_lo) in ions per cell or None
CAMPAIGNS = [
    ("JVASP-2017-LCO-B88", "LiCoO2, optB88-vdW+U", None),
    ("JVASP-2017-LCO-PBE", "LiCoO2, PBE+U", None),
    ("JVASP-2017-LCO", "LiCoO2, optPBE-vdW+U (original run)", None),
    ("JVASP-144791-NMC-vdW", "NMC cell, optB88-vdW+U", None),
    ("JVASP-144791-NMC", "NMC cell, PBE+U (original run)", None),
    # Lei et al. 2014: O3-Na1.00CoO2 charge curve over x = 1 -> 0.5
    ("JVASP-79809-NCO", "O3-NaCoO2, optB88-vdW+U", (8, 4)),
    # Okamoto et al. 2015 / Hatakeyama et al. 2019: Mg extraction up to x ~ 0.4
    ("JVASP-11340-MMO", "MgMn2O4, PBE+U", (16, 10)),
    ("JVASP-141543-B88", "Rb2LiFeF6 (A), optB88-vdW+U", None),
    ("JVASP-154749-B88", "Rb2LiFeF6 (B), optB88-vdW+U", None),
    ("JVASP-116849-B88", "LiV2F7, optB88-vdW+U", None),
    ("JVASP-77457-B88", "LiCa2Ag, optB88-vdW", None),
]


def hull_plateaus(pts, n_total, z, mu):
    """Lower convex hull over (n, E) points; returns plateaus as dicts."""
    pts = sorted(pts, key=lambda p: -p[0])           # from most to least ion
    e_hi, e_lo = pts[0][1], pts[-1][1]
    n_hi, n_lo = pts[0][0], pts[-1][0]
    span = n_hi - n_lo
    fe = [(n, e - ((n - n_lo) / span) * e_hi - ((n_hi - n) / span) * e_lo, e) for n, e in pts]
    hull = [fe[0]]
    for p in fe[1:]:
        while len(hull) > 1:
            (x0, f0, _), (x1, f1, _) = hull[-2], hull[-1]
            t = (x0 - x1) / (x0 - p[0])
            if f1 <= f0 + t * (p[1] - f0) + 1e-10:
                break
            hull.pop()
        hull.append(p)
    return [{"x_from": round(a[0] / n_total, 4), "x_to": round(b[0] / n_total, 4),
             "V": round((b[2] - a[2]) / (z * (a[0] - b[0])) + mu / z, 3)}
            for a, b in zip(hull, hull[1:])]


def converged(step_dir: Path):
    o = step_dir / "results" / "OUTCAR"
    if not o.exists():
        return None
    return "reached required accuracy" in o.read_text(errors="replace")


def main():
    out, rows = {}, []
    for name, label, rng in CAMPAIGNS:
        sup = next((ROOT / "dft_inputs" / name).glob("supercell_*"), None)
        if sup is None:
            continue
        data = json.load(open(sup / "energies.json"))
        ion, z, n_total, mu = dp._ion_info(data)
        if mu is None:
            mu = dp._E_ION_METAL.get(ion, {}).get(data.get("functional", "pbe"))
        E = {dp._step_n(s): s["energy"] for s in data["steps"] if s["energy"] is not None}
        if len(E) < 2 or mu is None:
            out[name] = {"label": label, "status": f"{len(E)} energies recorded" + ("" if mu is not None else ", no metal reference")}
            continue
        unconv = []
        for d in sup.glob("step_*"):
            m = dp._STEP_RE.match(d.name)
            if m and int(m.group(3)) in E and converged(d) is False:
                unconv.append(d.name)
        n_hi, n_lo = max(E), min(E)
        avg = lambda a, b: round((E[b] - E[a]) / (z * (a - b)) + mu / z, 3)
        rec = {"label": label, "ion": ion, "z": z, "functional": data.get("functional"),
               "mu_metal": mu, "n_total": n_total, "recorded": f"{n_hi}->{n_lo} of {n_total}",
               "V_avg_recorded_range": avg(n_hi, n_lo),
               "complete": n_hi == n_total and n_lo == 0,
               "hull": hull_plateaus(list(E.items()), n_total, z, mu),
               "unconverged_steps": sorted(unconv)}
        if rng and rng[0] in E and rng[1] in E:
            rec["V_avg_exp_range"] = {"range": f"{rng[0]}->{rng[1]}", "V": avg(*rng)}
        out[name] = rec
        try:  # plots only for tagged revision campaigns; never overwrite the published untagged figures
            if data.get("tag"):
                dp.compute_voltage_curve(str(sup), plot=True, e_ion_metal=mu)
        except Exception as e:  # plotting is secondary to the numbers
            rec["plot_error"] = str(e)
        rows.append(rec | {"name": name})

    json.dump(out, open(HERE / "revision_voltages.json", "w"), indent=1)
    lines = ["| Campaign | Functional | Range | V_avg (V) | Exp. range V_avg (V) | Hull plateaus (V) | Unconverged |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        er = r.get("V_avg_exp_range", {})
        lines.append(f"| {r['label']} | {r['functional']} | {r['recorded']} | {r['V_avg_recorded_range']:.2f} | "
                     f"{(er.get('range','') + ': ' + format(er['V'], '.2f')) if er else ''} | "
                     f"{', '.join(format(h['V'], '.2f') for h in r['hull'])} | {', '.join(r['unconverged_steps'])} |")
    (HERE / "revision_voltages.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for k, v in out.items():
        if "status" in v:
            print(f"{k}: {v['status']}")


if __name__ == "__main__":
    main()
