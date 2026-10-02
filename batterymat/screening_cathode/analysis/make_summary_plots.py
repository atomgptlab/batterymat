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
AVG_V_DFT    = [3.60, 3.91, 4.08, 4.40, 4.21]
MAX_V_ALIGNN = [3.74, 3.60, 4.42, 4.77, 5.06]
MAX_V_DFT    = [3.99, 4.34, 4.41, 5.03, 5.21]
CAP_ALIGNN   = [632, 606, 667, 702, 1372]
CAP_DFT      = [605, 586, 607, 679, 1380]  # LCO: relaxed optB88 step_00 volume 258.04 A^3

# Experimental references (see screening_cathode_analysis.md for citations).
# Avg: quasi-equilibrium average discharge voltage (GITT / low C-rate).
# Max: highest plateau voltage reported in galvanostatic cycling.
AVG_V_EXP = [3.43, 4.10, 4.05, 3.70, 4.05]  # LFP: Yamada 2001 two-phase equilibrium 3.43 V
MAX_V_EXP = [3.50, 4.10, 4.17, 4.20, 4.17]  # LFP: Padhi 1997 plateau 3.5 V

# Theoretical volumetric capacity (Q_grav [mAh/g] * density [g/cm^3]).
# LFP 170*3.60; LMP 171*3.43; LMO 148*4.28 (1 Li/f.u. reversible);
# NMC 278*4.77 (NMC-111 full delithiation); LCO 274*5.05.
CAP_EXP = [612, 587, 633, 1326, 1383]


def grouped_bar(ax, labels, groups, title, ylabel, ylim=None):
    names = list(groups.keys())
    values = [groups[n] for n in names]
    x = np.arange(len(labels))
    width = 0.8 / len(names)
    offsets = (np.arange(len(names)) - (len(names) - 1) / 2) * width
    colors = ["#4C78A8", "#F28E2B", "#59A14F", "#B07AA1"]
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
    ax.legend()


# Voltage summary (two panels)
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
grouped_bar(
    axes[0], MATERIALS,
    {"ALIGNN (tier 1)": AVG_V_TIER1, "ALIGNN-FF (tier 2)": AVG_V_ALIGNN, "DFT (tier 3)": AVG_V_DFT, "Experiment": AVG_V_EXP},
    "Average voltage",
    "Average Voltage (V)", ylim=(0, 6),
)
grouped_bar(
    axes[1], MATERIALS,
    {"ALIGNN-FF": MAX_V_ALIGNN, "DFT": MAX_V_DFT, "Experiment": MAX_V_EXP},
    "Maximum voltage",
    "Max Voltage (V)", ylim=(0, 8),
)
fig.tight_layout()
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
fig.tight_layout()
fig.savefig("capacity_summary.png", dpi=300)
plt.close(fig)

print("Wrote voltage_summary.png and capacity_summary.png")
