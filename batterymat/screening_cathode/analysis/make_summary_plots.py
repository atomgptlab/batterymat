"""Regenerate voltage_summary.png and capacity_summary.png with a third
Experimental column alongside ALIGNN-FF and DFT values.

Experimental voltages come from the literature compiled in
screening_cathode_analysis.md. Capacities use the theoretical crystallographic
volumetric capacity (Q_grav * density), consistent with the DFT/ALIGNN-FF bars
which are also theoretical (full delithiation on the relaxed/unrelaxed cell).
"""
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 14, "axes.labelsize": 15, "axes.titlesize": 15,
                     "legend.fontsize": 13, "xtick.labelsize": 14, "ytick.labelsize": 13})

MATERIALS = ["LFP", "LMP", "LMO", "NMC", "LCO"]

# Source provenance (regenerate these arrays when energies.json or Li_min.csv changes):
#   AVG_V_ALIGNN, MAX_V_ALIGNN: avg_voltage / max_voltage columns of
#     average_voltage/Li_min.csv for JVASP-{42723,116897,141792,144791,2017}.
#   AVG_V_DFT, MAX_V_DFT: mean / max of step voltages returned by
#     dft_prep.compute_voltage_curve on each material's energies.json.
#   CAP_ALIGNN: theoretical Q_grav * density on the unrelaxed JARVIS volume.
#   CAP_DFT: theoretical Q_grav * density on the relaxed CONTCAR volume.
AVG_V_ALIGNN = [3.49, 3.17, 4.07, 4.18, 3.84]
# Tier-1 ALIGNN scalar regressor (analysis/tier1/benchmarks_tier1.csv, v_tier1)
AVG_V_TIER1  = [3.25, 3.10, 3.60, 3.70, 3.42]
# LCO: optB88-vdW+U rerun (JVASP-2017-LCO-B88), average over x = 1 -> 0.5 (Li8 -> Li4), the
# experimentally cycled range; full-range average is 4.45 V (see revision_voltages.json).
AVG_V_DFT    = [3.60, 3.93, 4.08, 4.40, 4.21]
# LFP and LMP above are the JARVIS entries JVASP-42723/116897, a Pnma polymorph with tetrahedral Li.
# Olivine two-endpoint DFT averages with the protocol of the chains (Table 8 of the paper;
# dft_inputs/extra_2026-10-26/olivine-*/), drawn as markers on the DFT bars: index -> V.
AVG_V_DFT_OLIVINE = {0: 3.55, 1: 3.92}
MAX_V_ALIGNN = [3.74, 3.60, 4.42, 4.77, 5.06]
MAX_V_DFT    = [3.99, 4.34, 4.41, 5.03, 5.21]
CAP_ALIGNN   = [632, 606, 667, 702, 1372]
CAP_DFT      = [605, 586, 607, 679, 1380]  # LCO: relaxed optB88 step_00 volume 258.04 A^3

# Experimental references (order LFP, LMP, LMO, NMC, LCO).
# Avg: LFP 3.43 V two-phase equilibrium potential from GITT (Kondo et al., JES 165,
#   A2047, 2018); LMP ~4.1 V plateau (Padhi 1997; Aravindan et al., JMCA 1, 3518, 2013);
#   LMO 4.04 V, composition-weighted open-circuit voltage (Ohzuku, Kitagawa & Hirai,
#   JES 137, 769, 1990: 4.110 V for 0.27 < x < 0.60, taken to x = 0; one-phase region
#   centered at 3.94 V for 0.60 < x < 1); NMC 4.37 V, open-circuit voltage of NMC-111
#   averaged over Li0.5 -> Li0.03 per TM, the range the Li4Mn3Co2Ni3O16 cell covers
#   (Yabuuchi & Ohzuku, 206th ECS Meeting, Abstract 324, 2004; conventional-window
#   average 3.85 V, Yoshizawa & Ohzuku, Electrochemistry 71, 1177, 2003); LCO 3.97 V,
#   composition-weighted open-circuit voltage over x = 1 -> 0.5 (Ohzuku & Ueda, JES 141,
#   2972, 1994: 3.92 V plateau to x ~ 0.74, 4.20 V at x = 0.5; digitized curve).
# Max: highest plateau where the profile has one (LFP 3.5 V, Padhi 1997; LMP 4.1 V;
#   LMO 4.11 V, Ohzuku 1990); for the sloping NMC and LCO profiles the highest
#   open-circuit voltage over the compared range (NMC 4.71 V at Li0.03, Yabuuchi &
#   Ohzuku 2004; LCO 4.20 V at x = 0.5, Ohzuku & Ueda 1994).
AVG_V_EXP = [3.43, 4.10, 4.04, 4.37, 3.97]
MAX_V_EXP = [3.50, 4.10, 4.11, 4.71, 4.20]

# Theoretical volumetric capacity (Q_grav [mAh/g] * density [g/cm^3]).
# LFP 170*3.60; LMP 171*3.43; LMO 148*4.28 (1 Li/f.u. reversible);
# NMC 278*4.77 (NMC-111 full delithiation); LCO 274*5.05.
CAP_EXP = [612, 587, 633, 1326, 1383]


# one color per method in every panel
COLORS = {"ALIGNN (tier 1)": "#B07AA1", "ALIGNN-FF": "#4C78A8", "DFT": "#F28E2B", "Experiment": "#59A14F"}


def _color(name):
    return next(c for k, c in COLORS.items() if name == k or name.startswith(k + " "))


def grouped_bar(ax, labels, groups, title, ylabel, ylim=None, legend_kw=None):
    names = list(groups.keys())
    values = [groups[n] for n in names]
    x = np.arange(len(labels))
    width = 0.8 / len(names)
    offsets = (np.arange(len(names)) - (len(names) - 1) / 2) * width
    colors = [_color(n) for n in names]
    for off, name, vals, color in zip(offsets, names, values, colors):
        bars = ax.bar(x + off, vals, width, label=name, color=color)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + (0.04 if ylim else 8),
                    f"{v:.2f}" if ylim else f"{int(round(v))}",
                    ha="center", va="bottom", fontsize=10, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(**(legend_kw or {}))


# Voltage summary (two panels)
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
grouped_bar(
    axes[0], MATERIALS,
    {"ALIGNN (tier 1)": AVG_V_TIER1, "ALIGNN-FF (tier 2)": AVG_V_ALIGNN, "DFT (tier 3)": AVG_V_DFT, "Experiment": AVG_V_EXP},
    "Average voltage",
    "Average Voltage (V)", ylim=(0, 6.6), legend_kw={"loc": "upper center", "ncol": 2, "fontsize": 12},
)
grouped_bar(
    axes[1], MATERIALS,
    {"ALIGNN-FF": MAX_V_ALIGNN, "DFT": MAX_V_DFT, "Experiment": MAX_V_EXP},
    "Maximum voltage",
    "Max Voltage (V)", ylim=(0, 8),
)
# olivine markers on the DFT bars (third of four series, offset +0.1)
axes[0].scatter([i + 0.1 for i in AVG_V_DFT_OLIVINE], list(AVG_V_DFT_OLIVINE.values()), marker="D", s=55,
                facecolor="white", edgecolor="black", linewidth=1.3, zorder=5, label="DFT, olivine")
for t in axes[0].texts:  # lift the value labels of those bars clear of the markers
    if any(abs(t.get_position()[0] - (i + 0.1)) < 1e-6 for i in AVG_V_DFT_OLIVINE):
        t.set_y(t.get_position()[1] + 0.12)
axes[0].legend(loc="upper center", ncol=2, fontsize=12)
fig.text(0.005, 0.98, "a", fontsize=20, fontweight="bold", va="top")
fig.tight_layout(rect=(0.015, 0, 1, 1))
fig.savefig("voltage_summary.png", dpi=300)
plt.close(fig)

# Capacity summary (single panel)
fig, ax = plt.subplots(figsize=(10, 6))
grouped_bar(
    ax, MATERIALS,
    {"ALIGNN-FF": CAP_ALIGNN, "DFT": CAP_DFT, "Experiment": CAP_EXP},
    "Theoretical volumetric capacity",
    "Volumetric Capacity (mAh/cm³)",
)
ax.set_ylim(0, 1650)  # headroom for the rotated value labels of the LCO bars
fig.text(0.005, 0.98, "b", fontsize=20, fontweight="bold", va="top")
fig.tight_layout(rect=(0.02, 0, 1, 1))
fig.savefig("capacity_summary.png", dpi=300)
plt.close(fig)

print("Wrote voltage_summary.png and capacity_summary.png")
