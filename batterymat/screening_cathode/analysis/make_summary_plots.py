"""Regenerate voltage_summary.png and capacity_summary.png with a third
Experimental column alongside ALIGNN-FF and DFT values.

Experimental voltages come from the literature compiled in
screening_cathode_analysis.md. Capacities use the theoretical crystallographic
volumetric capacity (Q_grav * density), consistent with the DFT/ALIGNN-FF bars
which are also theoretical (full delithiation on the relaxed/unrelaxed cell).
"""
import matplotlib.pyplot as plt
import numpy as np

MATERIALS = ["LFP", "LMP", "LMO", "NMC", "LCO"]

# ALIGNN-FF and DFT values reproduced from the previous summary plots.
AVG_V_ALIGNN = [3.49, 3.17, 4.07, 4.18, 3.84]
AVG_V_DFT    = [3.60, 3.91, 4.08, 4.40, 4.18]
MAX_V_ALIGNN = [3.74, 3.60, 4.42, 4.77, 5.06]
MAX_V_DFT    = [3.99, 4.34, 4.41, 5.03, 4.78]
CAP_ALIGNN   = [632, 606, 667, 702, 1372]
CAP_DFT      = [605, 586, 607, 679, 1360]

# Experimental references (see screening_cathode_analysis.md for citations).
# Avg: quasi-equilibrium average discharge voltage (GITT / low C-rate).
# Max: highest plateau voltage reported in galvanostatic cycling.
AVG_V_EXP = [3.42, 4.10, 4.05, 3.70, 4.05]
MAX_V_EXP = [3.45, 4.10, 4.17, 4.20, 4.17]

# Theoretical volumetric capacity (Q_grav [mAh/g] * density [g/cm^3]).
# LFP 170*3.60; LMP 171*3.43; LMO 148*4.28 (1 Li/f.u. reversible);
# NMC 278*4.77 (NMC-111 full delithiation); LCO 274*5.05.
CAP_EXP = [612, 587, 633, 1326, 1383]


def grouped_bar(ax, labels, groups, title, ylabel, ylim=None):
    names = list(groups.keys())
    values = [groups[n] for n in names]
    x = np.arange(len(labels))
    width = 0.27
    offsets = np.linspace(-width, width, len(names))
    colors = ["#4C78A8", "#F28E2B", "#59A14F"]
    for off, name, vals, color in zip(offsets, names, values, colors):
        bars = ax.bar(x + off, vals, width, label=name, color=color)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}" if ylim else f"{int(round(v))}",
                    ha="center", va="bottom", fontsize=9)
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
    {"ALIGNN-FF": AVG_V_ALIGNN, "DFT": AVG_V_DFT, "Experiment": AVG_V_EXP},
    "Average Voltage: ALIGNN-FF vs DFT vs Experiment",
    "Average Voltage (V)", ylim=(0, 6),
)
grouped_bar(
    axes[1], MATERIALS,
    {"ALIGNN-FF": MAX_V_ALIGNN, "DFT": MAX_V_DFT, "Experiment": MAX_V_EXP},
    "Max Voltage: ALIGNN-FF vs DFT vs Experiment",
    "Max Voltage (V)", ylim=(0, 8),
)
fig.tight_layout()
fig.savefig("voltage_summary.png", dpi=150)
plt.close(fig)

# Capacity summary (single panel)
fig, ax = plt.subplots(figsize=(10, 6))
grouped_bar(
    ax, MATERIALS,
    {"ALIGNN-FF": CAP_ALIGNN, "DFT": CAP_DFT, "Experiment": CAP_EXP},
    "Volumetric Capacity: ALIGNN-FF vs DFT vs Experiment",
    "Volumetric Capacity (mAh/cm³)",
)
fig.tight_layout()
fig.savefig("capacity_summary.png", dpi=150)
plt.close(fig)

print("Wrote voltage_summary.png and capacity_summary.png")
