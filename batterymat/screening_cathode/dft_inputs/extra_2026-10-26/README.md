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
| `Li3FeO3-Li16-restart/` (added 2026-10-27) | `../JVASP-117295/supercell_1x1x2/step_08_Li16` | Li3FeO3, Li16Fe8O24 (48 atoms), end of the Fe3+/Fe4+ range | over the last 50 steps the energy spans 0.03 eV; no O-O bonds |

Each restart starts from the step-200 geometry (`results/CONTCAR` of the original step, used as
`POSCAR`) and from that run's final `WAVECAR` and `CHGCAR` (copied in on the cluster from
`/data/jlee859/dft/neurips/dft_inputs/<step>/`; not in the repo), with `ISTART = 1`,
`ICHARG = 1` as in the original INCAR. The INCAR is the as-run INCAR of the original step
with two changes: `IBRION = 2` (conjugate gradient, as in the earlier input template,
instead of the quasi-Newton `IBRION = 1` that oscillated) and `NSW = 300`. KPOINTS and
POTCAR_spec are unchanged (NELECT 256, 312 and 136, as in the original OUTCARs).

Intermediate steps that also hit NSW = 200 (Li3FeO3 steps 7 and 9-23, LiH4SNO4 steps 6-7) are not
restarted: they do not enter the average voltage.

Result of `Li3FeO3-Li0-restart` (2026-10-27): not usable. Continuing the relaxation of Fe8O24
formed four O2 molecules (O-O 1.23-1.25 A, so Fe8O24 -> Fe8O16 + 4 O2) and lowered the energy from
-158.59 to -174.38 eV by step 226, still falling; the conjugate-gradient search then diverged (step
231: O-O 0.67 A, forces up to 415 eV/A), so the final "reached required accuracy" and F = -159.31 eV
are spurious. The fully delithiated composition has no intercalation energy to report (the
step-200 energy of the original run, behind the earlier 4.26 V, predates the O2 formation). Li3FeO3
is therefore reported over the Fe3+/Fe4+ range, Li24 -> Li16 (steps 0 and 8), and step 8 is
restarted in `Li3FeO3-Li16-restart/` with the same two INCAR changes (its first energy should
reproduce the step-200 value, -260.70859 eV). Check its whole energy history and O-O distances, not
only the convergence line.

Result of `Li3FeO3-Li16-restart` (job 15357, 2026-10-27): converged in 5 conjugate-gradient steps,
energy falling monotonically by 8 meV to F = -260.71661191 eV (max force 0.09 eV/A), no O-O bond
(shortest 2.42 A), Fe moments 3.5 muB (Fe4+). Recorded as step 8 of JVASP-117295 (NSW-200 value kept as
`energy_nsw200`); Li3FeO3 over Li24 -> Li16: 3.82 V.


## Cell-relaxed endpoints of LFP, LMP and LMO

The LFP, LMP and LMO chains relax the cell only at full lithiation (step 0, ISIF = 3) and keep it
fixed afterwards (ISIF = 2). For the NMC variant, relaxing the cell at x = 0.5 and x = 0 changed the
averages by 0.21-0.36 V (`../extra_2026-10-01/NMC-{PBE,vdW}-ISIF3/`). These six runs repeat that
test for the other three PBE+U chains. The average over a composition range depends only on its two
endpoints, so with step 0 (already cell-relaxed) they give the cell-relaxed x = 1 -> 0.5, 0.5 -> 0
and 1 -> 0 averages directly.

| Directory | Starts from | Atoms |
|-----------|-------------|-------|
| `LFP-ISIF3/Li8`, `LFP-ISIF3/Li0` | `../JVASP-42723-LFP/supercell_2x2x1/step_08_Li8`, `step_16_Li0` | Li8Fe16P16O64, Fe16P16O64 |
| `LMP-ISIF3/Li8`, `LMP-ISIF3/Li0` | `../JVASP-116897-LMP/supercell_2x2x1/step_08_Li8`, `step_16_Li0` | Li8Mn16P16O64, Mn16P16O64 |
| `LMO-ISIF3/Li8`, `LMO-ISIF3/Li0` | `../JVASP-141792-LMO/supercell_2x2x2/step_08_Li8`, `step_16_Li0` | Li8Mn32O64, Mn32O64 |

Same protocol as the NMC test: the POSCAR is the final geometry of the chain step (rebuilt from its
OUTCAR, since the chain CONTCARs are not kept locally; the rebuild reproduces the LMP step-16 CONTCAR
to 5e-6 A and each step-8 geometry matches the step-9 POSCAR to 1e-5 A), and the INCAR is that
step's as-run INCAR with four changes: `ISIF = 3`, `NSW = 300`, and a cold start (`ISTART = 0`,
`ICHARG = 2`) so that a changing cell does not reuse a wave function from another cell. KPOINTS,
POTCAR_spec, MAGMOM and the U values are unchanged (NELECT 712, 688, 696, 672, 824, 800 as in the
chain OUTCARs). LMP and LMO step 8 used the earlier input template (IBRION = 2, NELM = 500, default
mixing), and it is kept. Check after the runs: the first ionic step is at the chain geometry and
cell, so its energy should reproduce the chain energy if the cold start lands in the same magnetic
state (chain energies: LFP -717.30743 / -672.96607, LMP -753.41416 / -704.89376, LMO -723.21011 /
-674.90141 eV for Li8 / Li0).

## Olivine endpoints of LiFePO4 and LiMnPO4

The LFP and LMP benchmark chains (JVASP-42723, mp-850917; JVASP-116897, mp-776296) are not the
olivine phase. Both are a Pnma polymorph with tetrahedral Li (Li-O 1.96-2.14 A; lattice 4.92 x
6.09-6.36 x 9.39 A; 24 and 15 meV/atom above the hull, no ICSD entry), whereas olivine has
octahedral Li (Li-O 2.10-2.21 A; 4.74 x 6.11 x 10.44 A for LiMnPO4). These four runs give the
average voltage of the olivine phases from their two endpoints, with the cell relaxed at both.

| Directory | Starts from | Atoms |
|-----------|-------------|-------|
| `olivine-LMP/Li16`, `olivine-LMP/Li0` | JVASP-91261 (olivine LiMnPO4, mp-18997, 1.5 meV/atom above the hull), 2x2x1; Li removed for Li0 | Li16Mn16P16O64, Mn16P16O64 |
| `olivine-LFP/Li16`, `olivine-LFP/Li0` | the same cell with Fe on the Mn sites (olivine LiFePO4 is not in the JARVIS snapshot) | Li16Fe16P16O64, Fe16P16O64 |

Built by `build_olivine.py` (reproducible from the JARVIS snapshot). The 2x2x1 supercell follows the
7 A rule (112 atoms, as in the benchmark chains). INCAR = `LFP-ISIF3/Li0/INCAR` (current template,
`ISIF = 3`, `NSW = 300`, cold start) with the LDAU tags and MAGMOM set per composition; U(Fe) 5.3 eV,
U(Mn) 3.9 eV; KPOINTS 2x2x3 as in the chains. Magnetic order: antiferromagnetic with zero net
moment, every nearest-neighbor pair of transition metals (corner-sharing octahedra in the same
layer, 3.92 A) antiparallel and every nearest pair in adjacent layers (5.51 A) antiparallel, the
same pattern in every primitive cell (k = 0). Expected NELECT: 736 / 688 (LFP Li16 / Li0), 720 /
672 (LMP Li16 / Li0). Average voltage over x = 1 -> 0: (E_Li0 - E_Li16)/16 + mu_Li with mu_Li =
-1.9031 eV/atom.

`olivine-LFP/Li0-fixcell` and `olivine-LMP/Li0-fixcell` (built by `build_olivine_fixcell.py`) hold the
delithiated endpoint at the relaxed Li16 cell (`ISIF = 2`, Li16 CONTCAR with all Li removed, otherwise
the Li0 inputs), so that the pair Li16 / Li0-fixcell follows the fixed-cell protocol of the benchmark
chains and the pair Li16 / Li0 gives the cell-relaxed value. The Li16 runs converged (30 and 15 ionic
steps), stayed Pnma with octahedral Li (Li-O 2.10-2.25 A) and kept the antiferromagnetic order
(Fe 3.76, Mn 4.65 muB); their lattice constants are about 1% above experiment, and they lie 29 (LFP) and
37 (LMP) meV/atom below the fully lithiated polymorph cells of the benchmark chains.

## Reruns of the three cancelled conjugate-gradient cell relaxations (2026-10-27)

`LiMnIr2-Li4-restart`, `LMP-ISIF3/Li8` and `LMO-ISIF3/Li8` (all `IBRION = 2` with `ISIF = 3`) stopped
writing output on 2026-10-26 between 18:22 and 18:47, each right after a converged electronic loop, and
were cancelled about seven hours later (jobs 15307, 15310, 15312). Before that, the electronic state of
LMP Li8 jumped at an unchanged geometry (steps 21-22 are 3.2 and 9.2 eV above step 20) and LMO Li8 reached
an unphysical energy 45 eV below its start in the electronic loop of step 11, which is the pattern of a
failed line search. Every `IBRION = 1` cell relaxation finished normally, so the reruns use it:

| Directory | Geometry | INCAR |
|-----------|----------|-------|
| `LMO-ISIF3/Li8-rerun` | ionic step 10 of `LMO-ISIF3/Li8` (F = -723.78786 eV, smooth up to there) | as `LMO-ISIF3/Li0` (current template: IBRION 1, NELM 100, AMIX/BMIX tags, ISIF 3, NSW 300, cold start) with the Li8 LDAU tags and MAGMOM |
| `LMP-ISIF3/Li8-rerun` | ionic step 20 of `LMP-ISIF3/Li8` (F = -755.23906 eV) | as `LMP-ISIF3/Li0` with the Li8 LDAU tags and MAGMOM |
| `LiMnIr2-Li4-rerun` | ionic step 30 of `LiMnIr2-Li4-restart` (F = -106.35547 eV) | as `LiMnIr2-Li4-restart` with IBRION 2 -> 1 and a cold start (ISTART 0, ICHARG 2) |

LiMnIr2: the original lithiated endpoint (IBRION 1, NSW 200) already fell from -99.10 to -101.75 eV, and the
restart fell a further 4.6 eV in 30 steps while the cell sheared (c 11.26 -> 9.68 A, b 11.27 -> 11.84 A,
volume 250 -> 238 A^3), so the structure is changing rather than finishing a relaxation. The delithiated
steps of the chain were relaxed in the old cell; once `LiMnIr2-Li4-rerun` converges, the delithiated
endpoint is recomputed at its cell (all Li removed, ISIF 2), as for the olivine pair, and the average
voltage is taken from that pair.

## Results of the cell-relaxed, olivine and rerun jobs (2026-10-27)

All runs below finished with "reached required accuracy" (checked over the whole energy history, not
only the last line), with the NELECT of the input and no O-O bond shorter than 2.40 A. Average
voltages use mu_Li = -1.9031 eV/atom and the cell-relaxed step 0 of each chain (LMP: the
alternating-order value).

Cell-relaxed endpoints (`LFP-ISIF3`, `LMP-ISIF3`, `LMO-ISIF3`; Li8 of LMP and LMO from the `-rerun`
directories):

| Chain | Li8: steps, F (eV) | Li0: steps, F (eV) | Volume vs Li16 (Li8 / Li0) | Fixed cell (chain), V | Cell relaxed, V |
|-------|--------------------|--------------------|----------------------------|------------------------|-----------------|
| LFP polymorph | 47, -719.00913991 | 37, -677.75143815 | -1.6% / -2.6% | 3.55 / 3.64 / 3.60 | 3.34 / 3.25 / 3.30 |
| LMP polymorph | 78, -755.27872898 | 222, -711.99260990 | -2.5% / -6.9% | 3.70 / 4.16 / 3.93 | 3.46 / 3.51 / 3.49 |
| LMO | 14, -723.73914922 | 24, -677.10991939 | -2.3% / -6.0% | 4.02 / 4.14 / 4.08 | 3.95 / 3.93 / 3.94 |

Voltages are x = 1 -> 0.5 / 0.5 -> 0 / 1 -> 0. The cold starts did not always land in the electronic
state of the chain at the same geometry: LFP Li8 started 0.405 eV above chain step 8, LMO Li8 0.103 eV
below it, and the Li8 reruns started 0.069 eV (LMO) and 0.041 eV (LMP) above the last good step of the
cancelled runs. These offsets move the half-range values by at most 0.05 V (LFP) and 0.01 V (LMO, LMP)
and do not enter the full-range values. Final states: LMO Li8 8 Mn3+ / 24 Mn4+ (3.86 / 3.06 muB),
total moment 0; LMP Li8 6 clear Mn3+ (3.88-3.95 muB), 8 Mn2+ (4.54-4.64) and two intermediate sites
(4.14, 4.31), antiferromagnetic.

Olivine endpoints (Li16 at ISIF 3: LFP 30 steps, F = -764.16915887 eV; LMP 15 steps, -802.41965953 eV):

| Material | Li0 cell relaxed (ISIF 3) | V | Li0 at the Li16 cell (ISIF 2) | V | Experiment |
|----------|---------------------------|---|-------------------------------|---|------------|
| LiFePO4 | 54 steps, -678.53622211, volume -4.8% | 3.45 | 26 steps, -676.88849074, P = -34 kB | 3.55 | 3.43 |
| LiMnPO4 | 53 steps, -711.01379980, volume -5.8% | 3.81 | 19 steps, -709.29623066, P = -18 kB | 3.92 | 4.10 |

All four stay Pnma with the antiferromagnetic order (Fe3+ 4.33, Mn3+ 3.88-3.91 muB). The fixed-cell
pair follows the protocol of the benchmark chains; the cell-relaxed pair is 0.10 V lower for both.

`LiMnIr2-Li4-rerun` (job 15360) did not converge: it used all 300 ionic steps and ended at
F = -106.01006636 eV with a maximum force of 0.16 eV/A. Its energy fell from -106.528 to -106.901 eV in
the first three steps, rose to -105.80 eV by step 6 while the forces kept falling (0.80 -> 0.50 eV/A
by step 10), and then wandered between -105.70 and -106.02 eV for the remaining 294 steps. The Mn
moments stayed at 4.3 muB (ferromagnetic) throughout, the Ir moments changed between +-0.17 muB, and
the total moment between 18.0 and 18.8 muB, so the energy surface the optimizer follows is not
continuous: the electronic state (Hubbard occupations on Mn, induced moments on Ir) changes between
ionic steps. No Li0 run is built from it.
