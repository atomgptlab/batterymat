"""Collect campaign-C spot-check DFT energies and compare with the tier-1 labels.

For every entry in ``dft_inputs/spotcheck/spotcheck_manifest.json`` this reads
the final total energy of ``lithiated/`` and ``delithiated/`` (from
``results/OUTCAR`` or ``results/OSZICAR``, falling back to the run directory
itself) and computes

    V_DFT = (E_delithiated - E_lithiated) / n_Li + E_Li_metal

with the same sign convention as ``dft_prep.compute_voltage_curve``. It prints
a DFT-vs-label table with MAE / ME / RMSE, writes ``spotcheck_results.csv``
next to this script, and reports "no results yet" when nothing has finished.

Usage (from screening_cathode/analysis/):  python spotcheck_collect.py
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SC_DIR = HERE.parent
MANIFEST = SC_DIR / "dft_inputs" / "spotcheck" / "spotcheck_manifest.json"
# Endpoints that hit NSW were restarted from their last structure (2026-10-01); a restart
# directory with results/OUTCAR replaces the original endpoint.
RESTART = SC_DIR / "dft_inputs" / "extra_2026-10-01" / "spotcheck-restart"


def endpoint_dir(jid: str, side: str) -> Path:
    r = RESTART / jid / side
    return r if (r / "results" / "OUTCAR").exists() else MANIFEST.parent / jid / side


def read_outcar_energy(path: Path):
    """Last 'free  energy   TOTEN' in OUTCAR (eV), or None."""
    energy = None
    with path.open(errors="replace") as fh:
        for line in fh:
            if "free  energy   TOTEN" in line:
                try:
                    energy = float(line.split()[-2])
                except (IndexError, ValueError):
                    pass
    return energy


def read_oszicar_energy(path: Path):
    """Last 'F=' (smeared free energy, matches OUTCAR TOTEN) in OSZICAR, or None."""
    energy = None
    pat = re.compile(r"F=\s*([-+0-9.Ee]+)")
    with path.open(errors="replace") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                try:
                    energy = float(m.group(1))
                except ValueError:
                    pass
    return energy


def endpoint_energy(run_dir: Path):
    """Energy from results/OUTCAR, results/OSZICAR, OUTCAR, OSZICAR (first found)."""
    for sub in ("results", "."):
        for name, reader in (("OUTCAR", read_outcar_energy), ("OSZICAR", read_oszicar_energy)):
            p = run_dir / sub / name
            if p.exists():
                e = reader(p)
                if e is not None:
                    return e, str(p.relative_to(run_dir))
    return None, None


def main():
    if not MANIFEST.exists():
        print(f"no manifest at {MANIFEST}; run spotcheck_prepare.py first")
        return 1
    man = json.loads(MANIFEST.read_text())
    rows = []
    for s in man["structures"]:
        base = MANIFEST.parent / s["jid"]
        e_l, src_l = endpoint_energy(endpoint_dir(s["jid"], "lithiated"))
        e_d, src_d = endpoint_energy(endpoint_dir(s["jid"], "delithiated"))
        row = dict(jid=s["jid"], formula=s["formula"], n_li=s["n_li"],
                   functional=s["functional"], e_li_metal=s["e_li_metal"],
                   label_V=s["label_voltage_V"], e_lith=e_l, e_delith=e_d,
                   dft_V=None, error_V=None, status="pending")
        if e_l is not None and e_d is not None:
            row["dft_V"] = (e_d - e_l) / s["n_li"] + s["e_li_metal"]
            row["error_V"] = row["dft_V"] - s["label_voltage_V"]
            row["status"] = f"done ({src_l}, {src_d})"
        elif e_l is not None or e_d is not None:
            row["status"] = "partial (one endpoint finished)"
        rows.append(row)

    done = [r for r in rows if r["dft_V"] is not None]
    n = len(rows)
    if not done:
        print(f"no results yet: 0/{n} spot-check pairs have final energies "
              f"(looking for results/OUTCAR or results/OSZICAR under {MANIFEST.parent})")
        partial = sum(1 for r in rows if r["status"].startswith("partial"))
        if partial:
            print(f"  {partial} structures have one endpoint finished")
        return 0

    print(f"{'JID':13s} {'formula':18s} {'func':9s} {'nLi':>3s} {'label V':>8s} {'DFT V':>8s} {'err V':>7s}")
    for r in rows:
        dv = f"{r['dft_V']:8.3f}" if r["dft_V"] is not None else "     ---"
        ev = f"{r['error_V']:+7.3f}" if r["error_V"] is not None else "    ---"
        print(f"{r['jid']:13s} {r['formula']:18s} {r['functional']:9s} {r['n_li']:3d} "
              f"{r['label_V']:8.3f} {dv} {ev}  {r['status']}")
    errs = [r["error_V"] for r in done]
    mae = sum(abs(e) for e in errs) / len(errs)
    me = sum(errs) / len(errs)
    rmse = (sum(e * e for e in errs) / len(errs)) ** 0.5
    print(f"\ncompleted {len(done)}/{n}:  MAE = {mae:.3f} V   ME = {me:+.3f} V   RMSE = {rmse:.3f} V")

    import csv
    out = HERE / "spotcheck_results.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
