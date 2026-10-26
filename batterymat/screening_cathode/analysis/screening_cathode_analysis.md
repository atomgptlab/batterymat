# Screening Cathode Analysis: DFT vs Experiment

This document compiles the two main DFT-vs-experiment comparisons for the five-material cathode screening campaign (LFP, LMP, LMO, NMC, LCO): (1) convex hull equilibrium voltages and (2) volumetric capacities, both with full literature references.

## Part 1: Equilibrium Hull Voltage Comparison

The convex hull (equilibrium) voltage represents the thermodynamic discharge profile — flat plateaus corresponding to two-phase coexistence regions between stable Li compositions. These are compared against experimental quasi-equilibrium measurements (low C-rate or GITT) from literature.

The hull is computed from formation energies: `dE(x) = E(x) - x*E(1) - (1-x)*E(0)`. The lower convex hull identifies thermodynamically stable compositions; between two hull vertices the system exists as a two-phase mixture at constant voltage. The number of hull plateaus equals the number of hull segments (vertices minus one).

---

## LFP (LiFePO₄, olivine, PBE+U)

**17/17 steps complete (Li₁₆→Li₀).**

| Plateau | x_Li range | DFT Hull | Experiment | Difference |
|---------|-----------|----------|------------|------------|
| 1 | 0.69 → 1.0 | 3.48 V | 3.43 V | +0.05 V |
| 2 | 0.56 → 0.69 | 3.58 V | 3.43 V | +0.15 V |
| 3 | 0.38 → 0.56 | 3.65 V | 3.43 V | +0.22 V |
| 4 | 0 → 0.38 | 3.67 V | 3.43 V | +0.24 V |
| **Weighted avg** | **0 → 1.0** | **3.60 V** | **3.43 V** | **+0.17 V** |

**Discussion:** Experimentally, LiFePO₄ exhibits a single flat plateau at ~3.4 V due to a first-order two-phase reaction between LiFePO₄ and FePO₄ with no stable intermediate solid solutions at room temperature. The DFT hull instead produces 4 plateaus (5 hull vertices), which is an artifact of the finite 2×2×1 supercell (16 Li atoms). The discrete Li-vacancy orderings in the supercell create artificial intermediate stable compositions that would not exist in the thermodynamic limit (infinite system or disordered solid solution).

The key observation is that all four plateaus are clustered within a narrow 0.19 V range (3.48–3.67 V), which is consistent with a single flat plateau broadened by finite-size effects. The systematic overestimate of +0.06 to +0.25 V relative to experiment is typical for PBE+U with U(Fe)=5.3 eV on olivine LiFePO₄.

The slight upward slope from 3.48 V (high x) to 3.67 V (low x) reflects stronger Fe-Fe interactions at low Li content in the finite supercell — in the thermodynamic limit, these would average out to the single ~3.4 V plateau.

**Experimental references:** Padhi, Nanjundaswamy & Goodenough, *J. Electrochem. Soc.* **144**, 1188 (1997): reversible Li extraction/insertion at 3.5 V. Yamada, Chung & Hinokuma, *J. Electrochem. Soc.* **148**, A224 (2001): flat two-phase plateau at 3.4 V. The 3.43 V reference is the two-phase equilibrium potential from GITT, Kondo, Sasaki, Barai & Srinivasan, *J. Electrochem. Soc.* **165**, A2047 (2018).

---

## LMP (LiMnPO₄, olivine, PBE+U)

**17/17 steps complete (Li₁₆→Li₀).** The fully delithiated endpoint (MnPO₄, Mn³⁺ high-spin d⁴) was rerun on 2026-10-26 (`dft_inputs/extra_2026-10-26/LMP-Li0/`): the earlier attempts had used the step-15 POTCAR (Li_sv, Mn_pv, P) on a Mn, P, O cell (NELECT 576 instead of 672), which is why that run sat ~190 eV above the expected energy. With the correct POTCAR it converged in 46 ionic steps (Mn ±3.97 μB, O ≤ 0.07 μB).

Step 0 converged in the block AFM order, steps 1–16 in the per-atom alternating order. Step 0 was recomputed in the alternating order (`extra_2026-10-26/LMP-Li16-altmag/`): 0.161 eV per cell lower. `energies.json` uses that value (the block-order energy is kept as `energy_block_order`).

| Plateau | x_Li range | DFT Hull |
|---------|-----------|----------|
| 1 | 0.75 → 1.0 | 3.63 V |
| 2 | 0.69 → 0.75 | 3.65 V |
| 3 | 0.56 → 0.69 | 3.79 V |
| 4 | 0.50 → 0.56 | 3.86 V |
| 5 | 0.44 → 0.50 | 3.98 V |
| 6 | 0.25 → 0.44 | 4.13 V |
| 7 | 0.19 → 0.25 | 4.21 V |
| 8 | 0 → 0.19 | 4.24 V |
| **Weighted avg** | **0 → 1.0** | **3.93 V** (exp ~4.1 V, −0.17 V) |

**Discussion:** the average underestimates the experimental ~4.1 V plateau by 0.17 V, the opposite sign from LFP (+0.17 V) under the same Materials Project U calibration; consistent with the oxide-fitted U being too small for Mn (and too large for Fe) in the phosphates. Like LFP, LMP is a two-phase system with a single flat plateau experimentally; the finite supercell splits it into several narrow hull plateaus. With the block-order step 0 the average would be 3.92 V.

**Experimental references:** ~4.1 V plateau: Padhi, Nanjundaswamy & Goodenough, *J. Electrochem. Soc.* **144**, 1188 (1997) (Mn³⁺/Mn²⁺ at 4.1 V); Aravindan, Gnanaraj, Lee & Madhavi, *J. Mater. Chem. A* **1**, 3518 (2013). Two-phase mechanism: Delacourt et al., *Chem. Mater.* **16**, 93 (2004); Li, Azuma & Tohda, *Electrochem. Solid-State Lett.* **5**, A135 (2002)

---

## LMO (LiMn₂O₄, spinel, PBE+U)

**17/17 steps complete (Li₁₆→Li₀).**

| Plateau | x_Li range | DFT Hull | Experiment (open circuit) | Difference |
|---------|-----------|----------|---------------------------|------------|
| Upper | 0 → 0.44 | 4.17 V | 4.11 V (two cubic phases, 0.27 < x < 0.60) | +0.06 V |
| Lower | 0.44 → 1.0 | 4.00 V | 3.94 V (midpoint of the one-phase region, 0.60 < x < 1) | +0.06 V |
| **Weighted avg** | **0 → 1.0** | **4.08 V** | **4.04 V** | **+0.04 V** |

**Discussion:** The hull reproduces the two-step 4 V profile of the spinel, with both plateaus about 0.06 V above the measured open-circuit values. The computed step lies at x = 0.44, the measured one near x = 0.60. The experimental average weights the open-circuit voltage by composition, with the 4.11 V plateau taken to extend from x = 0.27 (the lowest measured composition) to x = 0: (0.60 × 4.11 + 0.40 × 3.94) / 1.0 = 4.04 V.

**Experimental references:** Ohzuku, Kitagawa & Hirai, *J. Electrochem. Soc.* **137**, 769 (1990) (open-circuit voltage 5 h after each 3 h current step: "4.110(5) V for 0.27(2) < x < 0.60(3)", "S-shaped OCV curve [mid-point, 3.94(1) V] for 0.60(3) < x < 1.0"); Tomon et al., *Commun. Chem.* **5**, 54 (2022) (galvanostatic, "two distinct plateaus at around 4.0 and 4.15 V").

---

## NMC (Li₄Mn₃Co₂Ni₃O₁₆, layered variant, PBE+U)

**17/17 steps complete (Li₁₆→Li₀).** JARVIS-DFT stores JVASP-144791 as *P*1, so the spacegroup rule routes it to PBE+U; that chain is the benchmark value.

Per transition metal the cell runs from Li₀.₅TMO₂ (x = 1) to TMO₂ (x = 0), the upper half of the lithium range of NMC-111. The reference is therefore the measured open-circuit voltage of NMC-111 over that same half, not the average over the usual cycling window, which starts from full lithiation.

| Cell range | Li per TM | DFT (PBE+U, fixed cell) | Experiment (NMC-111 open circuit) | Difference |
|-----------|-----------|-------------------------|-----------------------------------|------------|
| x = 1 → 0.5 | 0.5 → 0.25 | 3.94 V | 4.19 V | −0.25 V |
| x = 0.5 → 0 | 0.25 → 0 (measured to 0.03) | 4.86 V | 4.58 V | +0.28 V |
| **x = 1 → 0** | **0.5 → 0 (measured to 0.03)** | **4.40 V** | **4.37 V** | **+0.03 V** |

Hull plateaus: 3.52 V (x 1 → 0.88), 3.97 V (0.88 → 0.81), 4.08 V (0.81 → 0.56), 4.18 V (0.56 → 0.50), 4.86 V (0.50 → 0).

**Discussion:** The averages agree, but only because the first half is 0.25 V low and the second half 0.28 V high: the computed rise between the halves (0.92 V) is more than twice the measured one (0.39 V). The excess rise appears in all four variants (`JVASP-144791-NMC-vdW/` and the cell-relaxed endpoints in `dft_inputs/extra_2026-10-01/`): optB88-vdW+U with fixed cell 4.28 / 5.24 / 4.76 V, PBE+U with relaxed endpoints 3.73 / 4.63 / 4.18 V, optB88-vdW+U with relaxed endpoints 4.02 / 4.88 / 4.45 V (x = 1 → 0.5 / 0.5 → 0 / 1 → 0). So it is not the dispersion treatment or the cell constraint. Whether it is a GGA+U artifact of the most oxidized states would need a hybrid-functional or meta-GGA calculation of the x = 1, 0.5 and 0 compositions.

Over the usual cycling window, from full lithiation (Li 1 → about 0.4 per TM), NMC-111 averages 3.85 V. The original submission compared the cell with ~3.7 V over that window, which the cell does not cover; that is where its +0.70 V residual came from.

**Experimental references:** Yabuuchi & Ohzuku, 206th Meeting of The Electrochemical Society, Abstract 324 (2004), https://www.electrochem.org/dl/ma/206/pdfs/0324.pdf. Its Fig. 3 gives open-circuit voltages from intermittent charge (1 h on, 5 h off) up to 4.71 V, reached at Li₀.₀₃ by charge count. The values above are integrated from the figure's data series, with Li content = 1 − Q / 277.85 mAh/g. Because composition is counted from charge, side reactions above about 4.5 V make the deepest compositions look more delithiated than they are. Yoshizawa & Ohzuku, *Electrochemistry* **71**, 1177 (2003) report "average discharge voltage of 3.85 V for the cell operation in 3.00–4.45 V" at C/10.

---

## LCO (LiCoO₂, layered, optB88-vdW+U)

**Reported run: `JVASP-2017-LCO-B88` (true optB88-vdW: GGA=BO, PARAM1=0.1833333333, PARAM2=0.22), 9/9 steps converged, no manual intervention.** The run in the original submission (`JVASP-2017-LCO`) used GGA=OR, which is optPBE-vdW, not optB88-vdW; it is kept for comparison only.

Experimentally LiCoO₂ is cycled only to about x = 0.5 (4.2 V charge cut-off), so the benchmark compares the DFT average over x = 1 → 0.5. The measured open-circuit voltage is flat at 3.92 V for 3/4 < x < 1 and at 4.50 V for 0 < x < 1/4 (Ohzuku & Ueda 1994). Over x = 1 → 0.5 it stays at 3.92 V down to x ≈ 0.74 and then rises to 4.20 V at x = 0.5. Weighted by composition (from the digitized curve, uncertainty about 0.02 V) it averages 3.97 V, which is the reference. The original submission used 4.05 V, the midpoint of 3.92 and 4.2 V, which is not a weighted average.

| Range | Experiment | optB88-vdW+U | PBE+U (same cell) | optPBE-vdW+U (original run) |
|-------|-----------|--------------|-------------------|-----------------------------|
| 3/4 < x < 1 (two-phase plateau) | 3.92 V | 4.18 V (+0.26) | 3.82 V (−0.10) | 4.18 V |
| x = 1 → 0.5 (cycled range) | 3.97 V | **4.21 V (+0.24)** | 3.81 V (−0.16) | 4.01 V |
| 0 < x < 1/4 (two-phase plateau) | 4.50 V | 4.78 V (+0.28) | 3.81 V (−0.69) | 4.47 V |
| Full range x = 1 → 0 | n/a | 4.45 V | 3.84 V | 4.18 V |

Hull plateaus (optB88): 4.18 / 4.23 / 4.62 / 4.78 V. Cell volume within 1% of the lithiated cell down to x = 3/8, then contracts gradually to −6.4% at CoO₂ (no blow-up).

**Discussion:** optB88-vdW+U reproduces the 0.58 V rise between the two experimental two-phase plateaus to within 0.02 V, with a nearly uniform offset of +0.26 to +0.28 V, and the +0.24 V residual over the cycled range matches that offset. PBE+U on the same supercell is flat (3.81–3.82 V) and misses the rise, and its cell expands by 3.5% on delithiation instead of contracting: without dispersion the interlayer binding of the Li-poor states is wrong, which is the case for routing layered hosts to the vdW functional. The functional gap on LCO is 0.62 V (full range) and 0.39 V (x = 1 → 0.5). The original optPBE run's 28% endpoint expansion did not recur with the correct functional. The uniform offset points to a systematic shift (U calibration or Li reference) rather than a deep-delithiation failure.

**Experimental references:** Ohzuku & Ueda, *J. Electrochem. Soc.* **141**, 2972 (1994) (open-circuit voltages 10 h after each 2 h current step: "constant at 3.92 V in the range between 0 and 70 mAh·g⁻¹", 4.20 V at Li₁/₂CoO₂; source of the composition-weighted 3.97 V); Tan et al., *Materials* **14**, 242 (2021) (4.2 V cut-off, ~half the theoretical capacity); Reimers & Dahn, *J. Electrochem. Soc.* **139**, 2091 (1992) (0.4 ≤ x ≤ 1 phase diagram).

---

## Cross-Material Summary

| Material | Structure | Functional | Hull plateaus | DFT avg V | Exp avg V | Error |
|----------|-----------|-----------|---------------|-----------|-----------|-------|
| LFP | Olivine | PBE+U | 4 (artifact of finite cell; should be 1) | 3.60 V | 3.43 V | +0.17 V |
| LMP | Olivine | PBE+U | 8 (finite-cell splitting of one plateau) | 3.93 V | ~4.1 V | −0.17 V |
| LMO | Spinel | PBE+U | 2 (measured profile has two steps) | 4.08 V | 4.04 V | +0.04 V |
| NMC | Layered | PBE+U | 5 | 4.40 V† | 4.37 V† | +0.03 V (−0.25 / +0.28 by half) |
| LCO | Layered | optB88-vdW | 4; reproduces the measured 0.58 V plateau rise | 4.21 V‡ | 3.97 V‡ | +0.24 V |

† Over Li₀.₅ → Li₀ per transition metal, the range the cell covers (experiment to Li₀.₀₃).
‡ Average over x = 1 → 0.5, the experimentally cycled range (full range 4.45 V).

### Key findings

1. **All five averages are within 0.25 V of experiment**, each compared over the lithium range its cell covers or over which it is cycled.

2. **LCO reproduces the shape of the voltage profile.** With the correct optB88-vdW functional, the 0.58 V rise between the two measured two-phase plateaus is reproduced to within 0.02 V at a nearly uniform +0.26 to +0.28 V offset. PBE+U is flat and misses the rise.

3. **LMO reproduces the two-step profile**, with both plateaus about 0.06 V above the measured open-circuit values. The computed step is at x = 0.44, the measured one near 0.60.

4. **LFP reproduces the flat plateau.** The four hull plateaus lie within 0.19 V of each other, consistent with a single flat plateau broadened by finite-size effects. The average is +0.17 V above the 3.43 V two-phase equilibrium potential (GITT, Kondo et al. 2018).

5. **LMP is complete.** The fully delithiated MnPO₄ endpoint (Mn³⁺, d⁴) converged once it was rerun with the correct POTCAR; the earlier failure was a setup error, not a PBE+U convergence problem. The full-range average, 3.93 V, underestimates experiment by 0.17 V.

6. **NMC agrees on average but not in shape.** Over the matching lithium range the first half is 0.25 V low and the second half 0.28 V high. The rise between the halves is more than twice the measured one, with both functionals and with or without cell relaxation.

### Systematic errors by functional

| Functional | Materials | Typical error | Direction |
|-----------|-----------|--------------|-----------|
| PBE+U | LFP, LMO | +0.04 to +0.17 V | Overestimates |
| PBE+U | LMP | −0.17 V | Underestimates (U too low for Mn²⁺/³⁺) |
| PBE+U | NMC (matched Li range) | −0.25 V (first half), +0.28 V (second half) | Rise between the halves too steep |
| optB88-vdW | LCO | +0.24 to +0.28 V | Uniform overestimate across the plateaus and the cycled range |
| PBE+U | LCO (comparison) | −0.10 to −0.69 V | Flat profile; misses the rise at low Li content |

---

## Part 2: Experimental Volumetric Capacity Reference

The volumetric capacity in `capacity_summary.png` compares ALIGNN-FF (unrelaxed JARVIS volumes), DFT (relaxed CONTCAR volumes), and the theoretical experimental value (Q_grav × ρ) as a third bar. Below are the experimentally reported values from the literature used for that column and for broader comparison.

**Formula:** Q_vol (Ah/L) = Q_grav (mAh/g) × ρ (g/cm³), where ρ is the crystallographic density of the fully lithiated phase.

| Material | Grav. cap. (mAh/g) | Density (g/cm³) | Theoretical Q_vol (Ah/L) | Practical Q_vol (Ah/L) | Practical cutoff |
|----------|------|------|------|------|------|
| LFP (LiFePO₄) | 170 | 3.60 | ~612 | ~510–590 | Full (x=0) |
| LMP (LiMnPO₄) | 171 | 3.43 | ~587 | ~250–430 | Kinetics-limited |
| LMO (LiMn₂O₄) | 148 | 4.28 | ~633 | ~450–550 | Cycled to x≈0 |
| NMC (NMC-111) | 278 (full) / 160 (practical) | 4.77 | ~1326 | ~760 | 4.3 V cutoff |
| LCO (LiCoO₂) | 274 (full) / 140 (practical) | 5.05 | ~1383 | ~710 | x≥0.5 (4.2 V) |

### Per-Material Notes

**LFP (LiFePO₄):** Theoretical 170 mAh/g corresponds to full Li⁺ extraction (Fe²⁺/Fe³⁺). Crystallographic density 3.60 g/cm³ (Pnma olivine) gives ~612 Ah/L theoretical. Practical electrode-level capacities reach ~510–590 Ah/L in carbon-coated nanoparticle formulations; bulk micron-scale LFP is kinetically limited. References: Padhi, Nanjundaswamy & Goodenough, *J. Electrochem. Soc.* **144**, 1188 (1997); Yamada, Chung & Hinokuma, *J. Electrochem. Soc.* **148**, A224 (2001).

**LMP (LiMnPO₄):** Theoretical 171 mAh/g, density 3.43 g/cm³ → ~587 Ah/L. Practical capacities are significantly lower (~250–430 Ah/L) due to poor electronic conductivity (~10⁻¹⁰ S/cm, orders of magnitude below LFP) and Jahn-Teller distortion in delithiated Mn³⁺O₆ octahedra. Requires heavy carbon coating and nanosizing. Reference: Delacourt et al., *Chem. Mater.* **16**, 93 (2004); Martha et al., *J. Electrochem. Soc.* **156**, A541 (2009).

**LMO (LiMn₂O₄):** Theoretical capacity 148 mAh/g (only 1 Li per formula unit cycles reversibly between 4 V plateaus; the 3 V plateau on further lithiation causes Jahn-Teller capacity fade). Density 4.28 g/cm³ → ~633 Ah/L theoretical. Practical ~450–550 Ah/L at room temp; significant capacity fade at elevated temperatures due to Mn dissolution. References: Thackeray, David, Bruce & Goodenough, *Mater. Res. Bull.* **18**, 461 (1983); Ohzuku, Kitagawa & Hirai, *J. Electrochem. Soc.* **137**, 769 (1990).

**NMC (Li[Ni,Mn,Co]O₂):** Capacity depends strongly on composition and cutoff voltage. NMC-111 theoretical 278 mAh/g (full delithiation), practical ~160 mAh/g at 4.3 V cutoff (x_Li ≈ 0.4 remaining). Density ~4.77 g/cm³ → ~760 Ah/L practical. Higher-Ni compositions (NMC-622, NMC-811) reach 180–220 mAh/g practical. The JARVIS structure (Li₄Mn₃Co₂Ni₃O₁₆, ≈ NMC-334) has no direct commercial equivalent. References: Ohzuku & Makimura, *Chem. Lett.* **30**, 642 (2001); Noh, Youn, Yoon & Sun, *J. Power Sources* **233**, 121 (2013); Manthiram, *Nat. Commun.* **11**, 1550 (2020).

**LCO (LiCoO₂):** Theoretical 274 mAh/g (full x=0 delithiation). Density 5.05 g/cm³ (highest among the five) → ~1383 Ah/L theoretical. In practice, commercial cells cycle only to x ≈ 0.5 (4.2 V cutoff) giving ~140 mAh/g, ~710 Ah/L — deeper delithiation causes H2→H3 phase transition and irreversible capacity loss. Recent high-voltage LCO (4.5 V) reaches ~180 mAh/g, ~910 Ah/L. References: Mizushima, Jones, Wiseman & Goodenough, *Mater. Res. Bull.* **15**, 783 (1980); Reimers & Dahn, *J. Electrochem. Soc.* **139**, 2091 (1992).

### Caveats

1. **Crystallographic vs electrode-level capacity.** All values above are crystallographic (single-crystal density × gravimetric capacity). Real electrodes include binder (PVDF, ~5–10 wt%), conductive carbon (~5–10 wt%), and porosity (~30–40%), which reduce the packed-electrode volumetric capacity by 2–3×.

2. **Theoretical vs practical.** Theoretical values assume full Li⁺ extraction (x=0), which is rarely achievable — either thermodynamically (LMP, LCO at high V) or kinetically (LFP, LMP at high C-rate). Practical values depend on voltage cutoff, temperature, and cycling rate.

3. **DFT voltages vs DFT capacities are independent comparisons.** The voltage accuracy (within ~0.2 V of experiment for 4 of 5 materials) does not necessarily translate to capacity accuracy — capacity depends on relaxed cell volume, which is more sensitive to functional choice (PBE overestimates volumes by ~2–5%, systematically underestimating volumetric capacity).
