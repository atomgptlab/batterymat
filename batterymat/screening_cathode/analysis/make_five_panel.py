"""Five-panel ALIGNN-FF vs DFT voltage-profile figure (main-text Fig.).

For each benchmark: DFT step voltages (from energies.json), the lower-convex-hull
equilibrium plateaus over the recorded lithium range, the ALIGNN-FF step-voltage
profile from the screening pass (Li_min.csv), and the open-circuit (hull-average)
voltage stated in the panel. Fonts >= 13 pt.

Usage: python make_five_panel.py   (writes five_panel_alignnff_vs_dft.png)
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 14,
                     "legend.fontsize": 11, "xtick.labelsize": 12, "ytick.labelsize": 12})

MATS = [("JVASP-42723", "LFP", "Li$_x$FePO$_4$", 3.43),
        ("JVASP-116897", "LMP", "Li$_x$MnPO$_4$", 4.10),
        ("JVASP-141792", "LMO", "Li$_x$Mn$_2$O$_4$", 4.05),
        ("JVASP-144791", "NMC", "Li$_x$Mn$_3$Co$_2$Ni$_3$O$_{16}$", 3.70),
        ("JVASP-2017", "LCO", "Li$_x$CoO$_2$", 4.05)]
# chain directory when it differs from <jid>-<abbrev>: LCO uses the optB88-vdW rerun
DIRS = {"LCO": "JVASP-2017-LCO-B88"}
# measured two-phase plateaus (x_lo, x_hi, V): Ohzuku & Ueda, J. Electrochem. Soc. 141, 2972 (1994)
EXP_PLATEAUS = {"LCO": [(0.75, 1.0, 3.92), (0.0, 0.25, 4.50)]}
# experimental average refers to this lithium range (LCO is cycled only to x ~ 0.5)
EXP_RANGE = {"LCO": 0.5}


def lower_hull(steps, e_li):
    """Hull plateaus over the recorded range (same construction as dft_prep)."""
    n_tot = max(s["n_li"] for s in steps)
    hi = [s for s in steps if s["n_li"] == n_tot][0]
    lo = min(steps, key=lambda s: s["n_li"])
    x_hi, x_lo = 1.0, lo["n_li"] / n_tot
    pts = []
    for s in steps:
        x = s["n_li"] / n_tot
        t = (x - x_lo) / (x_hi - x_lo)
        fe = s["energy"] - t * hi["energy"] - (1 - t) * lo["energy"]
        pts.append((x, s["n_li"], s["energy"], fe))
    pts.sort(key=lambda p: p[0], reverse=True)
    hull = [pts[0]]
    for p in pts[1:]:
        while len(hull) > 1:
            x0, _, _, f0 = hull[-2]
            x1, _, _, f1 = hull[-1]
            x2, _, _, f2 = p
            if x0 == x2:
                break
            fi = f0 + (x0 - x1) / (x0 - x2) * (f2 - f0)
            if f1 <= fi + 1e-10:
                break
            hull.pop()
        hull.append(p)
    plateaus = []
    for a, b in zip(hull[:-1], hull[1:]):
        dn = a[1] - b[1]
        plateaus.append((a[0], b[0], (b[2] - a[2]) / dn + e_li))
    ocv = sum((a - b) * v for a, b, v in plateaus) / sum(a - b for a, b, v in plateaus)
    return plateaus, ocv, x_lo


def main():
    li = pd.read_csv(os.path.join(ROOT, "average_voltage", "Li_min.csv"))
    li = li[li.name.str.startswith("Li_")]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.ravel()
    for ax, (jid, ab, host, vexp) in zip(axes, MATS):
        ej = glob.glob(os.path.join(HERE, "..", "dft_inputs", DIRS.get(ab, f"{jid}-{ab}"), "supercell_*", "energies.json"))[0]
        d = json.load(open(ej))
        steps = [s for s in d["steps"] if s["energy"] is not None]
        steps.sort(key=lambda s: -s["n_li"])
        n_tot = d["n_li_total"]
        e_li = d["e_li_metal"]
        xs, vs = [], []
        for a, b in zip(steps[:-1], steps[1:]):
            dn = a["n_li"] - b["n_li"]
            xs.append((a["n_li"] + b["n_li"]) / 2 / n_tot)
            vs.append((b["energy"] - a["energy"]) / dn + e_li)
        ax.plot(xs, vs, "o-", color="#1f77b4", lw=1.5, ms=5, label="DFT step voltage")
        plateaus, ocv, x_lo = lower_hull(steps, e_li)
        lab = "DFT equilibrium (hull)"
        for a, b, v in plateaus:
            ax.plot([b, a], [v, v], "--", color="#d62728", lw=2.2, label=lab)
            lab = None
        prof = [float(v) for v in li[li.name.str.contains("_" + jid + "_")].iloc[0].voltage_profile.split(";")]
        n = len(prof)
        xf = [1 - (i + 0.5) / n for i in range(n)]
        ax.plot(xf, prof, "s-", color="#ff7f0e", lw=1.2, ms=4, alpha=0.9, label="ALIGNN-FF step voltage")
        x_exp = EXP_RANGE.get(ab, 0.0)
        ax.plot([x_exp, 1.0], [vexp, vexp], color="0.35", ls=":", lw=1.5, label="Experiment (avg)")
        lab2 = "Experiment (two-phase plateau)"
        for x0, x1, v in EXP_PLATEAUS.get(ab, []):
            ax.plot([x0, x1], [v, v], color="0.2", lw=2.5, label=lab2)
            lab2 = None
        if ab in EXP_RANGE:
            E = {s["n_li"]: s["energy"] for s in steps}
            n_hi, n_lo = n_tot, round(n_tot * EXP_RANGE[ab])
            v_rng = (E[n_lo] - E[n_hi]) / (n_hi - n_lo) + e_li
        rng = "" if x_lo == 0 else f" ($x\\geq{x_lo:.2f}$)"
        ax.set_title(f"{ab}, {jid}")
        extra = f"\nDFT avg ($x\\geq{EXP_RANGE[ab]:.1f}$): {v_rng:.2f} V" if ab in EXP_RANGE else ""
        exp_lab = f"Experiment ($x\\geq{EXP_RANGE[ab]:.1f}$)" if ab in EXP_RANGE else "Experiment"
        ax.text(0.03, 0.04, f"OCV (DFT hull avg){rng}: {ocv:.2f} V{extra}\nALIGNN-FF avg: {sum(prof)/n:.2f} V\n{exp_lab}: {vexp:.2f} V",
                transform=ax.transAxes, fontsize=11, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
        ax.set_xlabel(f"$x$ in {host}")
        ax.set_ylabel("Voltage (V vs. Li/Li$^+$)")
        ax.set_xlim(-0.02, 1.02)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].axis("off")
    h, l = [], []
    for a in axes[:-1]:
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in l: h.append(hh); l.append(ll)
    axes[-1].legend(h, l, loc="center", fontsize=13, frameon=False)
    fig.tight_layout()
    out = os.path.join(HERE, "five_panel_alignnff_vs_dft.png")
    fig.savefig(out, dpi=300)
    print("wrote", out)


if __name__ == "__main__":
    main()
