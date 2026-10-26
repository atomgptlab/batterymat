# LMP follow-up runs (2026-10-26)

Two PBE+U runs that finish the LiMnPO4 chain in `../JVASP-116897-LMP/supercell_2x2x1/`.
Both use the step-15 INCAR (current template) with a cold start (`ISTART = 0`, `ICHARG = 2`)
so that MAGMOM sets the starting magnetic order.

| Directory | Structure | MAGMOM | Purpose |
|-----------|-----------|--------|---------|
| `LMP-Li0/` | step-15 relaxed positions (`step_15_Li1/CONTCAR`) with the last Li removed; Mn16P16O64, 96 atoms | alternating `4 -4 ...`, as in steps 1-15 | Fully delithiated endpoint (step 16). The earlier step-16 run used the step-15 POTCAR (Li_sv, Mn_pv, P) on a Mn, P, O cell (NELECT 576 instead of 672), so no endpoint energy exists. `POTCAR_spec` is now Mn_pv, P, O. |
| `LMP-Li16-altmag/` | step-0 relaxed geometry (`step_00_Li16/CONTCAR`) | alternating, as in steps 1-15 | Step 0 converged to the block order `4*4 4*-4 4*4 4*-4`; steps 1-15 converged to the per-atom alternating order. This run gives the Li16 energy in the same order as the rest of the chain. ISIF = 2 at the step-0 cell (the fixed cell of steps 1-15). |

Both runs use the fixed cell of steps 1-15, with lattice lines copied from `step_01_Li15/POSCAR`. NSW follows the template rule `min(300, max(200, 2 * n_atoms))`. Apart from NSW, ISTART, ICHARG, MAGMOM and (Li0 only) the LDAU species lists, every INCAR tag matches `step_15_Li1/INCAR`.

Run each with `run_single_generic.sh` (64 ranks, as the chain); outputs land in `results/`.
