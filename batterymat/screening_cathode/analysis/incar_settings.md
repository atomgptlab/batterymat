# INCAR Settings for DFT Cathode Calculations

Summary of VASP INCAR settings used for all five cathode materials in the sequential supercell delithiation campaign.

## Overview

| Material | JID | Supercell | Atoms | Functional | U_eff (eV) | KPOINTS | PAW (POTCAR) |
|----------|-----|-----------|-------|------------|------------|---------|--------------|
| LFP (LiFePO₄) | JVASP-42723 | 2×2×1 | 112 | PBE+U | Fe=5.3 | 2×2×3 | Li_sv, Fe_pv, P, O |
| LMP (LiMnPO₄) | JVASP-116897 | 2×2×1 | 112 | PBE+U | Mn=3.9 | 2×2×3 | Li_sv, Mn_pv, P, O |
| LMO (LiMn₂O₄) | JVASP-141792 | 2×2×2 | 112 | PBE+U | Mn=3.9 | 2×2×2 | Li_sv, Mn_pv, O |
| NMC (Li₄Mn₃Co₂Ni₃O₁₆) | JVASP-144791 | 2×2×1 | 112 | PBE+U | Mn=3.9, Co=3.32, Ni=6.2 | 2×2×3 | Li_sv, Mn_pv, Co, Ni_pv, O |
| LCO (LiCoO₂) | JVASP-2017 | 2×2×2 | 32 | optB88-vdW+U in the reported chain `JVASP-2017-LCO-B88` (the original `JVASP-2017-LCO` chain was optPBE-vdW+U, labelled optB88 until 2026-09-28) | Co=3.32 | 2×2×2 | Li_sv, Co, O |

## Where the settings live

Every INCAR under `dft_inputs/` is the file used for the run that produced the recorded
energy (checked against the INCAR echoed at the top of each OUTCAR). Read those files for
any single step; this note summarizes them.

## Two input templates

The template changed during the benchmark campaigns. Shared by both: `ENCUT = 520`,
`PREC = Accurate`, `ISMEAR = 0`, `SIGMA = 0.05`, `ISPIN = 2`, `LASPH = .TRUE.`,
`LMAXMIX = 4`, `ISYM = 0`, `LORBIT = 11`, `LREAL = Auto`, `NCORE = 8`,
`NSW = min(300, max(200, 2 * n_atoms))`.

| Tag | Earlier template | Current template |
|-----|------------------|------------------|
| IBRION | 2 (conjugate gradient) | 1 (quasi-Newton) |
| EDIFF | 1E-4 (1E-6 for LMP steps 0-2 and the original LCO steps 0-2) | 1E-4 |
| EDIFFG | default (10 x EDIFF); -0.03 on most of the 1E-6 steps | default (10 x EDIFF) |
| NELM / NELMIN | 500 / default | 100 / 4 |
| ISTART / ICHARG | cold start | 1 / 1 (restart from the previous step's WAVECAR/CHGCAR when present) |
| Mixing | default | AMIX = 0.2, BMIX = 0.0001, AMIX_MAG = 0.4, BMIX_MAG = 0.0001 |
| KPAR / NSIM | 2 / 8, except LMP steps 0-2 and the original LCO steps 0-2 (defaults) | 2 / 8 |

Per chain, the current template is used for 14 of 17 LFP steps, 14 of 17 NMC steps, 3 of 17
LMO steps, 3 of 17 LMP steps (steps 14-16; the step-16 endpoint was rerun with a cold start, see below),
all 9 steps of the reported LCO chain (`JVASP-2017-LCO-B88`) and every other campaign
(LCO-PBE, NMC-vdW, NaCoO2, MgMn2O4, prospective, spot-check). Exceptions:

- Four transitional steps (LFP and NMC step 2, LMO and LMP step 13) use IBRION = 1, NELM = 100
  and restarts but not the conservative mixing; they were continued after the template change.
- LFP step 1 uses the earlier template with reduced mixing (AMIX = 0.1, BMIX = 0.01,
  AMIX_MAG = 0.2, BMIX_MAG = 0.0001): its first run reached the 500-iteration SCF cap and
  stopped, and the step was continued from its last structure.
- The fully delithiated step of the original LCO chain (`JVASP-2017-LCO/step_08_Li0`) was
  rerun with ISIF = 2, a cold start, ALGO = All, TIME = 0.4, AMIN = 0.01, AMIX = 0.1,
  AMIX_MAG = 0.2, NELM = 200, NELMIN = 6, NSIM = 4 and MAGMOM 1/-1 on Co, after the ISIF = 3
  relaxation expanded the cell by about 28%.

## Step-Dependent Settings

| Setting | Step 0 (fully lithiated) | Steps 1+ (delithiated) |
|---------|-------------------------|------------------------|
| ISIF | 3 (full cell + ion relaxation) | 2 for hosts outside the layered spacegroup whitelist (LFP, LMP, LMO, and NMC, whose JARVIS entry is P1); 3 for whitelisted layered hosts (LCO) |
| NSW | `min(300, max(200, 2 * n_atoms))` | same rule |

**ISIF rationale:** Step 0 needs full cell relaxation to establish the equilibrium lattice. Steps 1+ use ISIF=2 to avoid Pulay stress artifacts from the incomplete basis set at changed volumes. Layered hosts keep ISIF=3 because the c-axis must relax to accommodate interlayer changes on delithiation; ISIF=2 can produce negative voltages in layered systems.

## DFT+U Parameters (Dudarev, LDAUTYPE=2)

U_eff values from the [Materials Project Hubbard U calibration](https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values), calibrated for oxide systems. LDAUL=2 (d-orbitals) for transition metals, LDAUL=-1 for all others. LDAUJ=0 for all (Dudarev formalism: only U-J matters).

### Per-Material LDAU Tags

**LFP** (elements: Li, Fe, P, O):
```
LDAUL = -1  2 -1 -1
LDAUU =  0  5.30  0  0
LDAUJ =  0  0  0  0
```

**LMP** (elements: Li, Mn, P, O):
```
LDAUL = -1  2 -1 -1
LDAUU =  0  3.90  0  0
LDAUJ =  0  0  0  0
```

**LMO** (elements: Li, Mn, O):
```
LDAUL = -1  2 -1
LDAUU =  0  3.90  0
LDAUJ =  0  0  0
```

**NMC** (elements: Li, Mn, Co, Ni, O):
```
LDAUL = -1  2  2  2 -1
LDAUU =  0  3.90  3.32  6.20  0
LDAUJ =  0  0  0  0  0
```

**LCO** (elements: Li, Co, O):
```
LDAUL = -1  2 -1
LDAUU =  0  3.32  0
LDAUJ =  0  0  0
```

## vdW-DF Tags (layered materials)

optB88-vdW (the JARVIS-DFT functional; `dft_prep.py` default for layered materials since 2026-09-28):

```
GGA      = BO             # optB88 exchange
PARAM1   = 0.1833333333
PARAM2   = 0.2200000000
LUSE_VDW = .TRUE.         # vdW-DF nonlocal correlation
AGGAC    = 0.0            # no gradient correction to correlation
```

The original LCO chain (`JVASP-2017-LCO/`, 4.18 V) was run with `GGA = OR` + `LUSE_VDW` +
`AGGAC = 0.0`, which is **optPBE-vdW**, not optB88-vdW (VASP wiki, "Nonlocal vdW-DF
functionals"). It is kept as `functional = "optpbevdw"`; the true optB88-vdW rerun is
`JVASP-2017-LCO-B88/`.

## MAGMOM (Antiferromagnetic Patterns)

Step 0 uses block-structured AFM patterns where consecutive groups of magnetic atoms alternate sign. Default magnetic moments per element: Mn=4.0 μB, Fe=5.0 μB, Co=0.6 μB, Ni=2.0 μB. Non-magnetic atoms (Li, O, P) get 0.

Steps 1+ inherit MAGMOM from the previous step's INCAR with the removed Li's entry dropped. This preserves the block structure rather than regenerating per-atom alternating signs (which caused SCF divergence in DFT+U). The LMP chain predates the inheritance: its steps 1-15 were started from a per-atom alternating pattern (`4 -4 4 -4 ...`) and converged in that order, whereas step 0 converged in the block order. `dft_inputs/extra_2026-10-26/LMP-Li16-altmag/` is step 0 recomputed in the alternating order (0.16 eV lower per cell).

### Step 0 MAGMOM Per Material

**LFP** (112 atoms: 16 Li, 16 Fe, 16 P, 64 O):
```
MAGMOM = 16*0  4*5 4*-5 4*5 4*-5  80*0
```

**LMP** (112 atoms: 16 Li, 16 Mn, 16 P, 64 O):
```
MAGMOM = 16*0  4*4 4*-4 4*4 4*-4  80*0
```

**LMO** (112 atoms: 16 Li, 32 Mn, 64 O):
```
MAGMOM = 16*0  8*4 8*-4 8*4 8*-4  64*0
```

**NMC** (112 atoms: 16 Li, 12 Mn, 8 Co, 12 Ni, 64 O):
```
MAGMOM = 16*0  4*4 4*-4 4*4  4*-0.6 4*0.6  4*-2 4*2 4*-2  64*0
```

**LCO** (32 atoms: 8 Li, 8 Co, 16 O):
```
MAGMOM = 8*0  8*0.6  16*0
```

Note: LCO's Co³⁺ is low-spin (t₂g⁶, S≈0), so the 0.6 μB moment is small and not strongly AFM-ordered.

## Li Metal Reference

Voltage calculations require the energy of metallic Li as a reference. Both use BCC Li (JVASP-913, Im-3m), Li_sv PAW, ENCUT=520, Gamma 17×17×17 k-mesh.

| Functional | E_li_metal (eV/atom) | Directory |
|-----------|---------------------|-----------|
| PBE | -1.9031 | `dft_inputs/JVASP-913-Li/Li_sv_PBE/` |
| optB88-vdW (GGA=BO, PARAM1/2) | -0.9778 | `dft_inputs/JVASP-913-Li/Li_sv/static/` (local); confirmation rerun in `Li_sv_optB88vdW/` |
| optPBE-vdW (GGA=OR) | -0.9646 | `dft_inputs/JVASP-913-Li/Li_sv_optPBEvdW/` |

**Note (corrected 2026-09-28):** an earlier version of this note called `GGA=BO` "optB86b" and
`GGA=OR` "optB88". That was backwards: `GGA=BO` with PARAM1 = 0.1833333333, PARAM2 = 0.22 is
optB88-vdW, `GGA=OR` is optPBE-vdW, and optB86b-vdW is `GGA=MK`. The `Li_sv/` run was the
correct optB88-vdW reference. The original LCO chain and the -0.9646 reference are a
consistent optPBE-vdW pair.

## LMP Step 16 (fully delithiated MnPO4)

The run that survives from the original chain used the step-15 POTCAR (Li_sv, Mn_pv, P) on a
Mn, P, O cell, so VASP read Mn as Li, P as Mn and O as P (NELECT 576 instead of 672). Its energy,
about 190 eV above the value expected from step 15, reflects that setup error, not a false
electronic minimum. The endpoint was rerun on 2026-10-26 (`dft_inputs/extra_2026-10-26/LMP-Li0/`, copied to
`JVASP-116897-LMP/supercell_2x2x1/step_16_Li0/`) from the step-15 relaxed structure with the
correct POTCAR (Mn_pv, P, O), the step-15 INCAR with a cold start (ISTART = 0, ICHARG = 2) and
the alternating MAGMOM of steps 1-15; it converged in 46 ionic steps at -704.89376 eV.

## TB-mBJ Static Settings (Optional)

For band gap and DOS calculations at selected delithiation steps. Generated via `python dft_prep.py static JVASP-XXXXX 0 1 2`.

```
IBRION = -1         # No ionic relaxation
NSW    = 0
ENCUT  = 520
EDIFF  = 1E-6
METAGGA = MBJ       # Tran-Blaha modified Becke-Johnson
LASPH  = .TRUE.
ICHARG = 2          # Atomic superposition
ALGO   = All
TIME   = 0.4
LORBIT = 11
NEDOS  = 2000       # Dense DOS mesh
NELM   = 1000       # More SCF steps (meta-GGA converges slowly)
AMIX   = 0.02       # Very gentle mixing
BMIX   = 0.001
AMIN   = 0.01
LWAVE  = .TRUE.
LCHARG = .TRUE.
ISPIN  = 2
LMAXMIX = 4
ISTART = 0
ISYM   = 0
NCORE  = 8
```

No DFT+U tags (MBJ is a separate functional). MAGMOM inherited from the relaxation step's INCAR.
