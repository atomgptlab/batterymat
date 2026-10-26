# Follow-up runs (2026-10-26)

## LMP

Two PBE+U runs that finish the LiMnPO4 chain in `../JVASP-116897-LMP/supercell_2x2x1/`.
Both use the step-15 INCAR (current template) with a cold start (`ISTART = 0`, `ICHARG = 2`)
so that MAGMOM sets the starting magnetic order.

| Directory | Structure | MAGMOM | Purpose |
|-----------|-----------|--------|---------|
| `LMP-Li0/` | step-15 relaxed positions (`step_15_Li1/CONTCAR`) with the last Li removed; Mn16P16O64, 96 atoms | alternating `4 -4 ...`, as in steps 1-15 | Fully delithiated endpoint (step 16). The earlier step-16 run used the step-15 POTCAR (Li_sv, Mn_pv, P) on a Mn, P, O cell (NELECT 576 instead of 672), so that run gives no usable endpoint energy. `POTCAR_spec` is now Mn_pv, P, O. |
| `LMP-Li16-altmag/` | step-0 relaxed geometry (`step_00_Li16/CONTCAR`) | alternating, as in steps 1-15 | Step 0 converged to the block order `4*4 4*-4 4*4 4*-4`; steps 1-15 converged to the per-atom alternating order. This run gives the Li16 energy in the same order as the rest of the chain. ISIF = 2 at the step-0 cell (the fixed cell of steps 1-15). |

Both runs use the fixed cell of steps 1-15, with lattice lines copied from `step_01_Li15/POSCAR`. NSW follows the template rule `min(300, max(200, 2 * n_atoms))`. Apart from NSW, ISTART, ICHARG, MAGMOM and (Li0 only) the LDAU species lists, every INCAR tag matches `step_15_Li1/INCAR`.

Run each with `run_single_generic.sh` (64 ranks, as the chain); outputs land in `results/`.

## Results (2026-10-26, atomgptlab jobs 15257 and 15258)

| Run | Ionic steps | Energy (eV) | Mn moments |
|-----|-------------|-------------|------------|
| `LMP-Li0` | 46, converged | -704.89375546 | +-3.97 (Mn3+), O <= 0.07 |
| `LMP-Li16-altmag` | 5, converged | -798.22301208 | +-4.64 (Mn2+) |

The alternating-order Li16 energy is 0.161 eV below the block-order step 0 (-798.06163257).
Both values are recorded in `../JVASP-116897-LMP/supercell_2x2x1/energies.json` (step 0 keeps
the block-order value as `energy_block_order`), and the Li0 inputs and outputs are also copied to
`../JVASP-116897-LMP/supercell_2x2x1/step_16_Li0/`. Full-range LMP average: 3.93 V (3.92 V with
the block-order step 0).

## Prospective endpoint restarts

Three endpoints of the prospective chains used all 200 ionic steps (NSW = 200) without
meeting the energy criterion (EDIFFG unset, so 10 * EDIFF = 1e-3 eV). The average voltage
of a chain depends only on its two endpoints, so these three are restarted to convergence.

| Directory | Restarts | Material | State at step 200 |
|-----------|----------|----------|-------------------|
| `Li3FeO3-Li0-restart/` | `../JVASP-117295/supercell_1x1x2/step_24_Li0` | Li3FeO3, delithiated endpoint Fe8O24 (32 atoms) | over the last 50 steps the energy jumps by up to 0.22 eV per step and spans 0.66 eV; max force 0.45 eV/A |
| `LiH4SNO4-Li0-restart/` | `../JVASP-97428/supercell_2x1x1/step_08_Li0` | LiH4SNO4, delithiated endpoint H32S8N8O32 (80 atoms) | over the last 50 steps the energy drifts down by 22 meV (at most 2 meV per step); max force 0.04 eV/A |
| `LiMnIr2-Li4-restart/` | `../JVASP-81829/supercell_1x2x2/step_00_Li4` | LiMnIr2 (control), lithiated endpoint Li4Mn4Ir8 (16 atoms, ISIF = 3) | over the last 50 steps the energy jumps by up to 0.08 eV per step and spans 0.09 eV; max force 0.2 eV/A, pressure about -2 kB |

Each restart starts from the step-200 geometry (`results/CONTCAR` of the original step, used as
`POSCAR`) and from that run's final `WAVECAR` and `CHGCAR` (copied in on the cluster from
`/data/jlee859/dft/neurips/dft_inputs/<step>/`; not in the repo), with `ISTART = 1`,
`ICHARG = 1` as in the original INCAR. The INCAR is the as-run INCAR of the original step
with two changes: `IBRION = 2` (conjugate gradient, as in the earlier input template,
instead of the quasi-Newton `IBRION = 1` that oscillated) and `NSW = 300`. KPOINTS and
POTCAR_spec are unchanged (NELECT 256, 312 and 136, as in the original OUTCARs).

Intermediate steps that also hit NSW = 200 (Li3FeO3 steps 7-23, LiH4SNO4 steps 6-7) are not
restarted: they do not enter the average voltage.

