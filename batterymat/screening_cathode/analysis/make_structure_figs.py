"""SI figure: the five benchmark supercells (step_00 POSCARs) as ball-and-stick projections."""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 13})
MATS = [("JVASP-42723-LFP", "LiFePO$_4$ (LFP), 2x2x1, 112 atoms"),
        ("JVASP-116897-LMP", "LiMnPO$_4$ (LMP), 2x2x1, 112 atoms"),
        ("JVASP-141792-LMO", "LiMn$_2$O$_4$ (LMO), 2x2x2, 112 atoms"),
        ("JVASP-144791-NMC", "Li$_4$Mn$_3$Co$_2$Ni$_3$O$_{16}$ (NMC variant), 2x2x1, 112 atoms"),
        ("JVASP-2017-LCO", "LiCoO$_2$ (LCO), 2x2x2, 32 atoms")]
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.ravel()
for ax, (d, title) in zip(axes, MATS):
    pos = sorted(glob.glob(os.path.join(HERE, "..", "dft_inputs", d, "supercell_*", "step_00_*", "POSCAR")))[0]
    at = read(pos, format="vasp")
    plot_atoms(at, ax, radii=0.45, rotation="10x,-20y,0z", show_unit_cell=2)
    ax.set_title(title, fontsize=12)
    ax.set_axis_off()
axes[-1].set_axis_off()
axes[-1].text(0.0, 0.5, "Colours follow the jmol convention:\nLi violet, O red, P orange,\nFe brown, Mn purple, Co pink, Ni green.\nStep-0 (fully lithiated) supercells\nas written by dft_prep init.", fontsize=13, va="center")
fig.tight_layout()
out = os.path.join(HERE, "benchmark_structures.png")
fig.savefig(out, dpi=300)
print("wrote", out)
