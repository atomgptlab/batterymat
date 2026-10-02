# Revision DFT campaigns (npj Computational Materials major revision, due 2026-10-06)

> **2026-09-28 correction.** Every "optB88-vdW" INCAR in this project used `GGA = OR`, which is
> optPBE-vdW. `dft_prep.py` now writes true optB88-vdW (`GGA = BO`, `PARAM1 = 0.1833333333`,
> `PARAM2 = 0.22`) for `functional = "optb88vdw"` and keeps the old tags as `"optpbevdw"`. All
> new vdW inputs below (NMC-vdW, NCO, `ref-*/*_optB88vdW`, 12 spot-check runs) were switched.
> New: `JVASP-2017-LCO-B88/` (true-optB88 rerun of the LCO chain; the original `JVASP-2017-LCO/`
> is now `functional = "optpbevdw"`), `JVASP-913-Li/Li_sv_optB88vdW/` (confirmation run of the
> -0.9778 eV/atom optB88-vdW Li reference in `Li_sv/static`; the GGA=OR run moved to
> `Li_sv_optPBEvdW/`). Runners set NCORE per job so it divides the ranks per k-group.

## Status and results (2026-10-02, all runs complete)

Runs on atomgptlab (`/data/jlee859/revision`, cpu partition, VASP 6 `vasp_std`). Results are
pulled back with the include-list rsync (energies.json, OUTCAR, CONTCAR, OSZICAR, INCAR, POSCAR,
KPOINTS, POTCAR_spec; never POTCAR) and summarised by `../analysis/revision_voltages.py`
(writes `analysis/revision_voltages.{json,md}`). Cluster-results cutoff: Fri 2026-10-03.

| Run | Status | Result |
|---|---|---|
| `JVASP-2017-LCO-B88` (true optB88-vdW+U) | complete, 9/9, no intervention; volume -6.4% at CoO2 | x 1->0.5: **4.21 V** (exp ~4.05, midpoint of 3.92 to 4.2 V cut-off); full 4.45 V; hull 4.18 / 4.23 / 4.62 / 4.78 V. Two-phase plateaus 4.18 / 4.78 vs exp 3.92 / 4.50 (Ohzuku & Ueda 1994): rise 0.59 vs 0.58 V (4.184 and 4.776 V unrounded), uniform +0.27 V offset |
| A. `JVASP-2017-LCO-PBE` | complete, 9/9; volume +3.5% at CoO2 relative to the relaxed LiCoO2 cell (+4.4% relative to the JARVIS input cell) | full 3.84 V, x 1->0.5 3.81 V, hull 3.73 / 3.80 / 3.98 V (flat, misses the rise). Functional gap 0.62 V full / 0.39 V matched |
| `JVASP-2017-LCO` (original, GGA=OR = optPBE-vdW) | kept for comparison | full 4.18 V, x 1->0.5 4.01 V |
| B. `JVASP-144791-NMC-vdW` | complete, 17/17 converged, no intervention (ISIF=2 after step 0, like the PBE chain) | full **4.76 V** (PBE+U 4.40; exp ~3.7); x 1->0.5 4.28 (PBE 3.94); x 0.5->0 5.24 (PBE 4.86). Near-uniform +0.36 V shift at fixed cell. Cell-relaxed endpoints (ISIF=3 at Li8 and Li0, `extra_2026-10-01/NMC-{PBE,vdW}-ISIF3/`): PBE+U 3.73 / 4.63 / 4.18 V (x 1->0.5 / 0.5->0 / full), optB88 4.02 / 4.88 / 4.45 V; cells contract 5% and 9-10%. So the fixed-cell protocol adds ~0.2 V, optB88 sits 0.25-0.38 V above PBE either way, and the ~0.9 V jump at x=0.5 persists in all four. PBE+U fixed-cell 4.40 V stays the benchmark (JARVIS spacegroup is P1, so the rule routes it to PBE) |
| C. `spotcheck/` | 28/30 converged (6 delithiated endpoints restarted from their last structure in `extra_2026-10-01/spotcheck-restart/`; Li8MnO6 and Li7FeO6 still drifting after 500 steps) | all: MAE 0.88 V, ME +0.68, Spearman 0.90; redox-metal hosts (n=14): MAE 0.54, offset +0.35, scatter 0.43, Spearman 0.94 (`analysis/spotcheck_stats.json`) |
| D. `JVASP-79809-NCO` (Na) | complete, 9/9 states; all converged except step 6 (Na2), which hit NSW = 200; cell -20% at CoO2 | Na8->Na4: **3.09 V** (exp ~2.9 over x 1->0.5, Lei et al. 2014), +0.19 V. Full range **3.73 V** (endpoints converged) against 3.17 V from the force-field screen (FF 0.56 V low); Na4->Na0 4.37 V; step voltages 3.32 / 3.00 / 3.13 / 2.89 / 4.82 / 4.17 / 3.55 / 4.95 V |
| D. `JVASP-11340-MMO` (Mg) | complete; steps 0 to 15 converged, step 16 (Mg0) hit NSW = 200 unconverged | Mg16->Mg10: **3.10 V** vs ~3.4 V (Okamoto et al. 2015), -0.30 V. Hull 3.08 / 3.11 / 3.32 / 3.43 / 3.66 V over Mg16->Mg1 (the unconverged Mg0 step excluded); Mg16->Mg1 3.22 V |
| D. metal references | complete | Na: PBE -1.31054, optB88 +0.91960; Mg: PBE -1.50576, optB88 +1.10509 eV/atom (in `dft_prep._E_ION_METAL`) |
| Li optB88 reference `JVASP-913-Li/Li_sv_optB88vdW` | complete | -0.97791 eV/atom (confirms -0.9778) |
| `JVASP-{141543,154749,77457,116849}-B88` (true-optB88 reruns of the prospective runs) | all complete | Rb2LiFeF6 A 6.127 V, B 6.128 V (optPBE 6.03); LiCa2Ag 0.446 V (optPBE 0.41); LiV2F7 4.918 V (optPBE 4.84). Prospective FF-vs-DFT MAE 0.88 -> 0.90 V |

All chains are complete and the local copies are current (2026-10-02).

VASP input bundles for four new campaigns answering the referee reports. Everything here
was generated on 2026-09-24 with `../dft_prep.py` (see "Changes to dft_prep.py" below) and
verified (POSCAR parses, atom and ion counts, MAGMOM length, POTCAR_spec coverage, KPOINTS,
functional tags, DFT+U values). No POTCAR files are included
(only `POTCAR_spec`); VASP outputs go to gitignored `results/` folders.

| Campaign | Referee item | Directory | Cell | Functional | VASP relaxations | Est. core-h (64 ranks) |
|---|---|---|---|---|---|---|
| A. LCO PBE+U series | R1.M4, R2.m3 | `JVASP-2017-LCO-PBE/supercell_2x2x2/` | 32 atoms, Li8 | PBE+U (U_Co 3.32) | 9 (Li8 -> Li0) | ~500 |
| B. NMC optB88-vdW+U series | R2.m5 | `JVASP-144791-NMC-vdW/supercell_2x2x1/` | 112 atoms, Li16 | optB88-vdW+U | 17 (Li16 -> Li0) | ~3,500 |
| C. Endpoint spot-check of ALIGNN labels | R1.M5, R2.M4 | `spotcheck/<JID>/{lithiated,delithiated}/` | primitive, 3 to 20 atoms | auto (PBE+U or optB88-vdW+U) | 60 (30 x 2) | ~1,500 to 3,000 |
| D. Na chain (O3-NaCoO2) | R1.M2 | `JVASP-79809-NCO/supercell_2x2x2/` | 32 atoms, Na8 | optB88-vdW+U | 9 (Na8 -> Na0) | ~500 |
| D. Mg chain (MgMn2O4 spinel) | R1.M2 | `JVASP-11340-MMO/supercell_2x2x2/` | 112 atoms, Mg16 | PBE+U (U_Mn 3.9) | 17 (Mg16 -> Mg0) | ~3,000 |
| D. metal references | R1.M2 | `ref-Na/Na_pv_{PBE,optB88vdW}/`, `ref-Mg/Mg_pv_{PBE,optB88vdW}/` | 1 to 2 atoms | both | 4 | < 10 |

Total: 116 relaxations, roughly 9,000 to 11,000 core-hours. Basis: the completed chains in this
folder ran at 64 MPI ranks; measured per-step means are 2.84 h (182 core-h) for the 112-atom
supercells (50 steps over LMP, LMO, NMC) and 0.98 h (63 core-h) for the 32-atom LCO
supercell (9 steps). optB88-vdW adds 10 to 30 % per step; the spot-check cells are small but
use denser k-meshes, so 30 to 50 core-h each is a reasonable planning number. Chains are
strictly serial per material (about 2 days wall time for a 32-atom chain, 2 to 3 weeks for a
112-atom chain at one step per 2 to 4 h), so launch all chains in parallel on day one.

## Directory tree (new entries only)

```
dft_inputs/
  README_revision_campaigns.md          this file
  JVASP-2017-LCO-PBE/supercell_2x2x2/   energies.json (tag LCO-PBE)  step_00_Li8/{POSCAR,INCAR,KPOINTS,POTCAR_spec}
  JVASP-144791-NMC-vdW/supercell_2x2x1/ energies.json (tag NMC-vdW)  step_00_Li16/...
  JVASP-79809-NCO/supercell_2x2x2/      energies.json (ion Na, z 1)  step_00_Na8/...
  JVASP-11340-MMO/supercell_2x2x2/      energies.json (ion Mg, z 2)  step_00_Mg16/...
  ref-Na/Na_pv_PBE/  ref-Na/Na_pv_optB88vdW/    bcc Na, JVASP-14608
  ref-Mg/Mg_pv_PBE/  ref-Mg/Mg_pv_optB88vdW/    hcp Mg, JVASP-14840
  spotcheck/spotcheck_manifest.json
  spotcheck/JVASP-XXXXXX/lithiated/    spotcheck/JVASP-XXXXXX/delithiated/   (30 JIDs)
../analysis/spotcheck_selection.csv   selection with rationale (campaign C)
../analysis/spotcheck_prepare.py      regenerates the selection and inputs (fixed seed 20260924)
../analysis/spotcheck_collect.py      DFT-vs-label table + MAE once results/ exist
../chain_step_generic.sh              SLURM chain driver (copy of chain_step.sh with a config block)
../run_single_generic.sh              SLURM single-directory runner (spot-check, references)
```

Existing directories (`JVASP-2017-LCO`, `JVASP-144791-NMC`, the other benchmark chains,
`JVASP-913-Li`) were not modified.

## Launching on a generic SLURM cluster

1. Copy to the cluster: `dft_prep.py`, `chain_step_generic.sh`, `run_single_generic.sh`,
   and `dft_inputs/` (only the new campaign folders are needed) into one work directory,
   e.g. `~/dft/revision/`. `dft_prep.py next` needs a conda env with `jarvis-tools`,
   `alignn` and `torch` (CPU is fine; the ALIGNN weights `alignnff_wt01` must be present or
   downloadable, see CLAUDE.md for the figshare WAF workaround).
2. Edit the `EDIT ME` block at the top of both scripts (WORKDIR, conda, VASP binary or
   module, POTCAR library, optional vdw_kernel.bindat) and the `#SBATCH` header (partition,
   account, ntasks, walltime). The scripts also accept the same settings as `BMAT_*`
   environment variables. Rank counts are scaled to the ion count and kept a multiple of
   KPAR (2), as in `chain_step.sh`.
3. Chains (A, B, D): pass the directory name, not the bare JID.
   ```
   sbatch --job-name=bm_JVASP-2017-LCO-PBE   chain_step_generic.sh JVASP-2017-LCO-PBE
   sbatch --job-name=bm_JVASP-144791-NMC-vdW chain_step_generic.sh JVASP-144791-NMC-vdW
   sbatch --job-name=bm_JVASP-79809-NCO      chain_step_generic.sh JVASP-79809-NCO
   sbatch --job-name=bm_JVASP-11340-MMO      chain_step_generic.sh JVASP-11340-MMO
   ```
   Each segment runs one step, records the energy, generates the next step (ALIGNN vacancy
   ranking), and resubmits itself; completed steps are skipped, so a failed segment is
   resumed by resubmitting the same command. `check_chains.sh` works unchanged for status.
4. Single runs (C and the references): one job per directory.
   ```
   for d in dft_inputs/spotcheck/*/*/; do n=$(echo $d | awk -F/ '{print $3"_"$4}'); sbatch --job-name="sc_$n" run_single_generic.sh "$d"; done
   for d in dft_inputs/ref-*/*/;        do sbatch --job-name="ref_$(basename $d)" run_single_generic.sh "$d"; done
   ```
   Use 8 to 16 ranks for these small cells (the runner caps ranks at the ion count).
5. Bring back to the repo: every `energies.json`, every `results/{OUTCAR,CONTCAR,OSZICAR}`
   (gitignored, kept locally), and the generated `step_XX` input folders. Then
   `python dft_prep.py voltage <dir name>` for the chains and
   `python analysis/spotcheck_collect.py` for campaign C.

## Campaign A: LCO with PBE+U (same 2x2x2 supercell as the optB88-vdW chain)

`python dft_prep.py init JVASP-2017 --functional pbe --suffix LCO-PBE --db-cache <dft_3d.pkl>`

`step_00_Li8/POSCAR`, `KPOINTS` (2 2 2) and `POTCAR_spec` (Li_sv, Co, O) are byte-identical to
`JVASP-2017-LCO/supercell_2x2x2/step_00_Li8/`; MAGMOM (`8*0 8*0.6 16*0`) and the LDAU block
(U_Co 3.32) are identical. ISIF=3 at step 0 and at every later step (layered, spacegroup 166,
`layered: true` in energies.json), exactly like the optB88-vdW chain, so the only physical
difference is the functional. Voltages use `e_li_metal = -1.9031` (PBE Li_sv reference).

INCAR differences versus the existing chain other than the removed vdW tags are template
drift: the old chain was generated with an earlier `_PBE_INCAR` (IBRION=2, EDIFF=1E-6,
EDIFFG=-0.03, NELM=500, ISTART=0, no KPAR/NSIM/mixing tags); the current template writes
IBRION=1, EDIFF=1E-4 (EDIFFG default), NELM=100, NELMIN=4, ISTART=1/ICHARG=1, KPAR=2, NSIM=8,
AMIX=0.2, BMIX=1E-4, AMIX_MAG=0.4, BMIX_MAG=1E-4. None of these change the converged
energy beyond a few meV per cell. For a strictly like-for-like comparison, copy the old
chain's INCAR tags (minus GGA/LUSE_VDW/AGGAC) over the new INCAR before launching; the
chain then propagates them (`next` regenerates INCARs from the current template, so apply
the same edit to each new step or accept the drift, which is what the NMC PBE chain and the
LCO vdW chain already differ by).

## Campaign B: NMC with optB88-vdW+U (same 2x2x1 supercell as the PBE chain)

`python dft_prep.py init JVASP-144791 --functional optb88vdw --suffix NMC-vdW --db-cache <dft_3d.pkl>`

POSCAR, KPOINTS (2 2 3), POTCAR_spec (Li_sv, Mn_pv, Co, Ni_pv, O), MAGMOM
(`16*0 4*4 4*-4 4*4 4*-0.6 4*0.6 4*-2 4*2 4*-2 64*0`) and the LDAU block (Mn 3.9, Co 3.32,
Ni 6.2) are identical to `JVASP-144791-NMC/supercell_2x2x1/step_00_Li16/`. Added tags:
`GGA = BO`, `PARAM1 = 0.1833333333`, `PARAM2 = 0.2200000000`, `LUSE_VDW = .TRUE.`, `AGGAC = 0.0`
(true optB88-vdW; switched from `GGA = OR` on 2026-09-28). The remaining INCAR differences are the same
template drift listed under campaign A. JVASP-144791 is spacegroup P1 in JARVIS, so
`layered` is false and later steps use ISIF=2 (like the PBE chain); if you want the layered
treatment (ISIF=3 throughout, as for LCO) set `"layered": true` in `energies.json` before the
first `next`. Voltages use `e_li_metal = -0.9778` (optB88-vdW Li_sv reference). Note the
stoichiometry caveat for this JID (Li:TM = 1:2, not NMC-111).

## Campaign C: DFT endpoint spot-check of the tier-1 labels (30 structures)

Selection (`analysis/spotcheck_prepare.py`, seed 20260924, saved with rationale in
`analysis/spotcheck_selection.csv`):

1. tier-1 held-out test ids (`average_voltage/Li_250_ids_train_val_test.json`, key `id_test`,
   761 ids; identical copy in the Li_250 archive),
2. intersected with the Li screening pool (`Li_min.csv` rows named `Li_*`, 7,193): 738,
3. primitive cell <= 20 atoms, ehull <= 0.10 eV/atom, no f-block elements, label in 0 to 6 V: 296,
4. stratified 5 per 1-V bin (bin counts 47/41/66/106/30/6, so no fill was needed).

Label = `avg_voltage_V` from the archive `Li_250/summary.csv` (the tier-1 training label).
Note that this label differs from the `avg_voltage` in `Li_min.csv` (older screen.py run)
by up to 1.7 V for some entries; both are stored in the selection CSV.

Inputs: primitive JARVIS cell, no supercell. `lithiated/` is the JARVIS structure,
`delithiated/` has every Li removed from the same cell. Both are ISIF=3 relaxations with the
`dft_prep.py` INCAR template, functional auto-selected as `dft_prep` does (optB88-vdW for
spacegroups 166/194/12/15, PBE otherwise; 6 of 30 are vdW), DFT+U and AFM MAGMOM from
`_afm_magmom()`. KPOINTS are Gamma-centred with n_i = ceil(21 A / |a_i|), which reproduces the
3x3x3 mesh `dft_prep` uses at its 7 A supercell threshold (meshes are 3 to 6 per axis here).
`spotcheck/spotcheck_manifest.json` lists jid, formula, n_li, functional, e_li_metal, label.

Average voltage from the two endpoints, same sign convention as `compute_voltage_curve`:

    V_DFT = (E_delithiated - E_lithiated) / n_Li + E_Li_metal     (E_Li_metal: -1.9031 PBE, -0.9778 optB88-vdW)

`python analysis/spotcheck_collect.py` reads `results/OUTCAR` (last TOTEN) or
`results/OSZICAR` (last F=) under each endpoint, prints the DFT-vs-label table with MAE, ME
and RMSE, writes `analysis/spotcheck_results.csv`, and prints "no results yet" until then.

Caveats to state in the response: the endpoint average ignores intermediate phases (the
label is also a rigid-lattice endpoint quantity, so the comparison is like-for-like); the
fully delithiated cells of alkali-rich compounds (e.g. Li8MnO6 -> MnO6, Li7FeO6 -> FeO6) may
relax far from the host framework under ISIF=3, and intermetallics/high-voltage fluorides
were kept because they are in the screening pool, not because they are realistic cathodes.

## Campaign D: Na and Mg chains

Chosen from the JARVIS dft_3d snapshot (criteria: contains the ion, ehull <= 0.05 eV/atom,
supercell <= 120 atoms under the 7 A rule, experimentally characterised chemistry):

| Ion | JID | Formula | Spacegroup | ehull (eV/atom) | Supercell | Why |
|---|---|---|---|---|---|---|
| Na | JVASP-79809 | NaCoO2 | 166 (R-3m), O3 | 0.000 | 2x2x2, 32 atoms, 8 Na | Direct Na analogue of the LCO chain (same R-3m cell, same U_Co, auto optB88-vdW, ISIF=3); O3-NaCoO2 is the stoichiometric x = 1.00 phase of Lei et al., Chem. Mater. 2014 (10.1021/cm5021788), cycled 2.5 to 3.4 V for 1.00 > x > 0.52. JVASP-1858 is the same structure at ehull 0.0006. |
| Mg | JVASP-11340 | MgMn2O4 | 141 (I4_1/amd), tetragonal spinel | 0.017 | 2x2x2, 112 atoms, 16 Mg | The Jahn-Teller distorted (hausmannite-type) spinel that is the experimental MgMn2O4 phase (Kim et al., Adv. Mater. 2015, 10.1002/adma.201500083; Okamoto et al., Adv. Sci. 2015, 10.1002/advs.201500072, Mg redox couples near 3.4 V (Mn3+/Mn4+) and 2.3 V (Mn2+/Mn3+) vs Mg). The cubic Fd-3m entry JVASP-11520 has ehull 0.058 (fails the cutoff); JVASP-10644 (R3m, ehull 0.0) and JVASP-10051 (Imma, 0.010) are not the spinel. |

Alternatives that also pass the filters: O'3-NaMnO2 JVASP-1861 (C2/m, ehull 0.009, 3x3x2 =
72 atoms, 18 Na); olivine NaFePO4 is not in JARVIS (only Cmcm maricite-like entries
JVASP-143994/48383). Experimental average voltages for the manuscript must be verified against
the cited papers before use; the values above are what the searches on 2026-09-24 returned.

Generation:
```
python dft_prep.py init JVASP-79809 --ion Na --suffix NCO --db-cache <dft_3d.pkl>   # optB88-vdW+U, ISIF=3 all steps
python dft_prep.py init JVASP-11340 --ion Mg --suffix MMO --db-cache <dft_3d.pkl>   # PBE+U, ISIF=2 after step 0
```
MAGMOM: NCO `8*0 8*0.6 16*0`; MMO `16*0 8*4 8*-4 8*4 8*-4 64*0` (AFM Mn blocks mapped from the
primitive cell, as for Li). POTCAR_spec uses the JARVIS labels Na_pv, Mg_pv. Capacities printed
by `init` (all ions removed, z = 1 for Na, z = 2 for Mg): NaCoO2 235 mAh/g, MgMn2O4 271 mAh/g.

Metal references. `ref-Na/Na_pv_*` (bcc Na, JVASP-14608, 1 atom, k 17x17x17) and
`ref-Mg/Mg_pv_*` (hcp Mg, JVASP-14840, 2 atoms, k 17x17x11) mirror `JVASP-913-Li/Li_sv_*`:
ENCUT 520, ISMEAR=1, SIGMA=0.1, ISPIN=1, EDIFF=1E-6, LREAL=.FALSE.; the optB88-vdW folders are
copies of `Li_sv_optB88vdW/INCAR` (ISIF=3 relaxation, EDIFFG=-1E-3). The PBE folders are the
same ISIF=3 relaxation without the vdW tags, rather than the static run used for
`Li_sv_PBE`, because that static used an MP PBE-relaxed geometry that is not available for
Na/Mg (the JARVIS cells are optB88-vdW-relaxed). E_ion_metal = final TOTEN / number of atoms
in the cell (divide by 2 for hcp Mg). Then

```
python dft_prep.py voltage JVASP-79809-NCO --e-ion-metal <E_Na, optB88-vdW>
python dft_prep.py voltage JVASP-11340-MMO --e-ion-metal <E_Mg, PBE>
```
or set `BMAT_E_ION_METAL` for the chain script; until then `energies.json` carries
`"e_ion_metal": null` and `voltage` raises a clear error. Once known, also fill
`_E_ION_METAL["Na"/"Mg"]` in `dft_prep.py`. Voltage with z electrons per ion:
`V = (E_lo - E_hi) / (z * dn) + E_metal / z`.

## Changes to dft_prep.py (2026-09-24)

All Li behaviour is byte-identical: `init` for JVASP-2017, JVASP-141792 and JVASP-144791 and
`next`/`record`/`voltage` on an LCO step were run with the pristine and the modified module
and `diff -r` of the output trees is empty (`scratchpad/hpc/checks/byte_identical.sh`).
New POSCAR/KPOINTS/POTCAR_spec also match the committed step_00 files; INCARs differ from the
committed ones only by the template drift described under campaign A, which predates this work.

- `--ion {Li,Na,K,Mg,Ca,Zn,Al}` on `init` (default Li). `next`, `record`, `static`, `voltage`
  read the ion from `energies.json`. Step directories are `step_XX_<Ion>N`. Vacancy ranking
  (`get_next_vacancy(atoms, ion=)`), MAGMOM inheritance and capacity use the ion. Non-Li
  `energies.json` carries `ion`, `z`, `n_ion_total`, `e_ion_metal`, and `n_ion` /
  `removed_ion_index` per step; Li files keep the legacy keys and both are read by
  `_ion_info()` / `_step_n()`.
- Voltage with z electrons per ion in `compute_voltage_curve` (step and hull voltages);
  `--e-ion-metal` on `init` and `voltage` (alias of `--e-li-metal`). `theoretical_capacity()`
  (gravimetric and volumetric, printed by `init`).
- `--suffix` on `init` writes `<JID>-<suffix>/` and stores `tag` in `energies.json`; plots are
  then named `voltage_curve_<JID>-<tag>.png` so the LCO-PBE and NMC-vdW curves do not overwrite
  `voltage_curve_JVASP-2017.png` / `voltage_curve_JVASP-144791.png`.
- `--db-cache <pkl|json>` on `init` (or env `BATTERYMAT_DFT3D_CACHE`) loads a local dft_3d
  snapshot instead of the figshare download.
- Directory lookup (`_find_supercell_dir`): an exact directory-name match now wins; prefix
  matching on the bare JID is unchanged when it is unique and raises a clear error otherwise
  (`JVASP-2017` is ambiguous between `JVASP-2017-LCO` and `JVASP-2017-LCO-PBE`; pass the full
  name). `chain_step.sh` already passes directory names, so it is unaffected.
- `compute_dft_voltage()` (unused legacy helper) had the metal-reference sign opposite to
  `compute_voltage_curve`; it now uses the same convention and takes `z`.
