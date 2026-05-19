# INCAR Settings for DFT Cathode Calculations

Summary of VASP INCAR settings used for all five cathode materials in the sequential supercell delithiation campaign.

## Overview

| Material | JID | Supercell | Atoms | Functional | U_eff (eV) | KPOINTS | PAW (POTCAR) |
|----------|-----|-----------|-------|------------|------------|---------|--------------|
| LFP (LiFePO₄) | JVASP-42723 | 2×2×1 | 112 | PBE+U | Fe=5.3 | 2×2×3 | Li_sv, Fe_pv, P, O |
| LMP (LiMnPO₄) | JVASP-116897 | 2×2×1 | 112 | PBE+U | Mn=3.9 | 2×2×3 | Li_sv, Mn_pv, P, O |
| LMO (LiMn₂O₄) | JVASP-141792 | 2×2×2 | 112 | PBE+U | Mn=3.9 | 2×2×2 | Li_sv, Mn_pv, O |
| NMC (Li₄Mn₃Co₂Ni₃O₁₆) | JVASP-144791 | 2×2×1 | 112 | PBE+U | Mn=3.9, Co=3.32, Ni=6.2 | 2×2×3 | Li_sv, Mn_pv, Co, Ni_pv, O |
| LCO (LiCoO₂) | JVASP-2017 | 2×2×2 | 32 | optB88-vdW+U | Co=3.32 | 2×2×2 | Li_sv, Co, O |

## Common Settings (All Materials)

These tags are shared across all five materials and all delithiation steps:

```
ENCUT  = 520        # Plane-wave cutoff (eV)
EDIFF  = 1E-6       # Electronic convergence (eV)
EDIFFG = -0.03      # Ionic convergence (eV/Å, force-based)
PREC   = Accurate
ISMEAR = 0          # Gaussian smearing
SIGMA  = 0.05       # Smearing width (eV)
ISPIN  = 2          # Spin-polarized
LASPH  = .TRUE.     # Non-spherical PAW contributions
LMAXMIX = 4         # l-mixing for d-electrons
ISYM   = 0          # No symmetry (broken by vacancies)
LORBIT = 11         # Projected DOS
LREAL  = Auto       # Real-space PAW projection (safe for supercells >20 atoms)
LWAVE  = .TRUE.     # Write WAVECAR
LCHARG = .TRUE.     # Write CHGCAR
IBRION = 2          # Conjugate gradient ionic relaxation
ISTART = 0          # Fresh wavefunction (no WAVECAR reuse)
```

**Note:** Step_00 INCARs use `ISTART=0` (cold start). Later steps can use `ISTART=1` / `ICHARG=1` to reuse WAVECAR/CHGCAR from the previous step, but the INCARs on disk have `ISTART=0` throughout.

## Step-Dependent Settings

| Setting | Step 0 (fully lithiated) | Steps 1+ (delithiated) |
|---------|-------------------------|------------------------|
| ISIF | 3 (full cell + ion relaxation) | 2 (ion-only, fixed cell shape) |
| NSW | 200--224 (material-dependent) | Decreases by 2 per step (fewer atoms) |

**ISIF rationale:** Step 0 needs full cell relaxation to establish the equilibrium lattice. Steps 1+ use ISIF=2 to avoid Pulay stress artifacts from the incomplete basis set at changed volumes. For layered materials (LCO), ISIF=3 should ideally be kept at all steps because the c-axis must relax to accommodate interlayer expansion on delithiation -- ISIF=2 can produce negative voltages in layered systems.

## Performance Tags

```
NCORE = 8           # Band parallelization (sqrt of 64 cores)
KPAR  = 2           # k-point parallelization (2 groups)
NSIM  = 8           # Bands optimized simultaneously
NELM  = 500         # Max SCF iterations
```

KPAR=2 and NSIM=8 are present in LFP, LMO, NMC, and LCO but absent from LMP (no effect on results, only performance). Adjust NCORE and KPAR for different core counts.

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

## optB88-vdW Tags (LCO Only)

Appended to the PBE base INCAR for layered LiCoO₂:

```
GGA      = OR       # optB88 exchange functional
LUSE_VDW = .TRUE.   # Enable vdW-DF nonlocal correlation
AGGAC    = 0.0      # Turn off GGA correlation (replaced by vdW-DF)
```

## MAGMOM (Antiferromagnetic Patterns)

Step 0 uses block-structured AFM patterns where consecutive groups of magnetic atoms alternate sign. Default magnetic moments per element: Mn=4.0 μB, Fe=5.0 μB, Co=0.6 μB, Ni=2.0 μB. Non-magnetic atoms (Li, O, P) get 0.

Steps 1+ inherit MAGMOM from the previous step's INCAR with the removed Li's entry dropped. This preserves the block structure rather than regenerating per-atom alternating signs (which caused SCF divergence in DFT+U).

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
| optB88-vdW (GGA=OR) | -0.9646 | `dft_inputs/JVASP-913-Li/Li_sv_optB88vdW/` |

**Note:** The original optB88-vdW Li reference (`dft_inputs/JVASP-913-Li/Li_sv/`, now removed) was computed with `GGA=BO` (optB86b exchange, -0.9778 eV/atom), while LCO cathode steps use `GGA=OR` (optB88 exchange). A corrected reference with `GGA=OR` gives -0.9646 eV/atom — a 0.013 eV difference. All LCO voltages now use the corrected value.

## LMP Step 16 (Abandoned)

The fully delithiated MnPO₄ (step_16, Li₀) required special cold-start settings due to Mn⁴⁺ (d³) convergence pathology. Multiple attempts failed -- the energy was ~200 eV above the expected ground state, indicating a trapped false electronic minimum from near-degenerate magnetic configurations.

Final INCAR attempt:
```
IBRION = 1          # RMM-DIIS (changed from CG)
NSW    = 300
ISIF   = 2
EDIFF  = 1E-4       # Loosened from 1E-6
NELM   = 100        # Fail-fast on bad SCF
NELMIN = 4
ISTART = 0          # Cold start (no WAVECAR)
ICHARG = 2          # Atomic superposition (respects MAGMOM)
ALGO   = All        # Davidson + RMM-DIIS hybrid
AMIX     = 0.1      # Gentler charge mixing
BMIX     = 3E-5
AMIX_MAG = 0.2      # Gentler magnetic mixing
BMIX_MAG = 3E-5
MAGMOM = 3 -3 3 -3 3 -3 3 -3 3 -3 3 -3 3 -3 3 -3  80*0
```

This step was abandoned as not experimentally relevant (olivine LMP is never fully delithiated in practice).

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
