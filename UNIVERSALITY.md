# qEW or qKPZ below the crossover length l₀?

Summary of `step10_universality.py` over all eleven studies, 2026-10-13; Studies 8 and 9 added 2026-10-14, Studies 10 (GdFeCo) and 11 2026-10-15. Definitions and options are in the [README](README.md) (step 10).

## Question
In d = 2 a driven domain wall has several regimes of scale. At the smallest scales the regime is non-universal; the Larkin regime would follow, but here it is below the optical resolution. Between that and the crossover length l₀, the creep avalanches are depinning-like. Above l₀, overhangs change the geometry and ζ_eff ≈ 0.5 (Kolton, Ferrero & Rosso, [arXiv:2306.13415](https://arxiv.org/abs/2306.13415); l₀ from `step8`).

The question is which depinning class describes the intermediate range:
- **qEW:** ζ ≈ 1.25, and λ = 0 or λ ∝ v.
- **qKPZ:** ζ ≈ 0.63, and λ stays finite as v → 0.

Two features of the data matter:
- The creep equilibrium regime (ζ = 2/3) is numerically close to qKPZ.
- The optical PSF inflates ζ at ℓ ≲ 4 µm (step9 self-test).

Two independent tests are therefore used, both measured on the walls themselves.

## Tests
**A. Roughness exponent in [4 µm, l₀).**
- Measured on the local width with local rotation (step7), folds excluded; the contour and u_area are cross-checks.
- Free power law over [4 µm, l₀].
- Broken power law of step8 over [4 µm, L] (ζ₂ = 0.5, l₀ free), with ζ₁ = 1.25, 0.63, 2/3 or free; χ² per model.
- Errors: block bootstrap over time.

**B. KPZ term.**
- Normal velocity of contour segments of arc length L = 2–16 µm (< l₀) vs their local orientation. The switched pixels of the arrival map are assigned to the nearest contour point, and segments near defects are dropped.
- Λ = λ/v is the curvature of v(m) = v₀[1 + (Λ/2) m²], the velocity projected on the mean direction of motion:

  | Λ | Meaning |
  |---|---|
  | 0 | elastic line driven along h (qEW) |
  | 1 | isotropic medium; the geometric term λ = v vanishes with v |
  | ≫ 1 | qKPZ: λ/v diverges in creep, where v → 0 |

- R(L) = |Λ| ⟨m_L²⟩/2 is the relative change of velocity produced by the tilts of scale L.
- Two estimates: the contour (TIFF studies) and the segments of h(x, t) (all studies).

## Results

| study | sample, drive | v [nm/s] | l₀ [µm] (step8, local width) | ζ, [4 µm, l₀], folds excluded | ζ₁, broken power law | χ² qEW / qKPZ (folds excl.) | Λ contour: median (range) | Λ h(x): median (range) | max R |
|---|---|---|---|---|---|---|---|---|---|
| 5 | Pt/Co/Pt (thesis), 0.36 A, 51 °C | 0.13 | 80 (64–80) | 1.15 (1.13–1.17) | **1.26 (1.24–1.28)** | 1108 / 10651 | 1.25 (0.54–1.37) | 1.21 (1.14–1.51) | 0.23 |
| 6 | same, 0.38 A, 51 °C | 0.12 | 13 (5–13) | 0.93 (0.91–0.97) | 0.76 (no crossover, poor fit) | 968 / 882 | 0.24 (−1.79–1.10) | −0.03 (−0.46–0.01) | 0.12 |
| 7 | same, 0.42 A, 52 °C | 0.57 | 15 (8–15) | 1.04 (1.02–1.06) | 1.18 (1.14–1.28) | 32 / 675 | 0.70 (0.47–1.02) | 0.27 (0.13–0.40) | 0.12 |
| 2 | same, 0.8 A, RT | 1.65 | 28 (8–28) | 1.06 (1.01–1.10) | 1.10 (1.04–1.14) | 34 / 1054 | 0.83 (0.64–1.20) | 0.57 (0.33–0.72) | 0.17 |
| 4 | Pt/Co/Pt (2016), 0.50 A | 0.43 | 15 (5–18) | 1.11 (1.06–1.15) | 1.40 (1.21–1.56) | 8 / 176 | 0.80 (0.23–1.40) | 0.95 (0.14–1.26) | 0.13 |
| 3 | same, 0.48 A | 0.21 | 10 (7–40) | 0.92 (0.88–0.95) | 1.09 (1.02–1.16) | 17 / 161 | 0.10 (−0.33–0.81) | −0.12 (−0.91–0.25) | 0.06 |
| 1 | video (other sample region), 25 fps | — | 28 (14–28) | 1.05 (1.03–1.07) | 1.05 (1.03–1.08) | 27 / 756 | — | −0.05 (−0.21–0.28) | 0.01 |
| 8 | P7925, 0.04 A (scale assumed) | 0.30 | 23 (14–23) | 1.06 (1.03–1.08) | 1.04 (0.99–1.11) | 412 / 2740 | 1.41 (0.26–1.66) | 0.83 (−0.45–1.04) | 0.25 |
| 9 | same, two days later | 0.19 | 48 (29–48, poorly bounded) | 0.90 (0.89–0.91) | ≈ 1.0 (no break in the field) | 1161 / 3948 | 1.22 (0.73–1.99) | 0.60 (0.48–0.75) | 0.37 |
| 10 | GdFeCo B1743, 600 mA (scale assumed) | 2.36 | 35 (18–35) | 1.17 (1.14–1.21) | 1.31 (1.26–1.37) | 13 / 1266 | 0.56 (0.10–2.57) | −0.02 (−1.48–0.36) | 0.24 |
| 11 | P7925, 0.04 A, wall barely moves (scale assumed) | 0.15 | 26 (5–26, no clean break) | 0.77 (0.76–0.79) | 0.89 (no break in the field) | 1700 / 3714 | 0.65 (0.07–1.23) | 0.09 (−0.44–0.61) | 0.17 |

Notes on the table:
- **Errors:** 16–84 % block-bootstrap intervals.
- **l₀:** the main step8 value, with the range of its three estimators in brackets.
- **Λ ranges:** only the (L, dt) values whose 16–84 % interval is narrower than 1.5. The excluded ones are a few points with ±2–8 error bars:
  - L = 16 µm in Study 5;
  - L = 8 µm in Study 6;
  - L = 8 µm with the longest lag in Study 3.
- **χ² values:** with a diagonal error model; neighbouring scales are correlated, so the Δχ² overstate the significance. The ordering is robust.

Figures: `fig_universality.png` in each study; Studies 1–7 in `0_42A_20s_20x_52C_1/fig_universality_compare.png`, Studies 8–9 in `P7925_0_04A_24jun16/fig_universality_compare.png`.

## Reading
1. **Roughness.** ζ ≈ 1.0–1.15 below l₀ with a free power law, and ζ₁ ≈ 1.1–1.4 when the crossover to ≈ 0.5 is fitted. No study is compatible with 0.63 or with 2/3.
   - The self-test of the local width (README, step7) shows that it underestimates ζ near 1, so the values are compatible with qEW's 1.25.
   - **Study 5 is the cleanest case:** its l₀ is large (≈ 70 µm), and the wall is self-affine with ζ = 1.25 over more than a decade.
   - **Weak cases:**
     - Study 6: l₀ ≈ 5–13 µm leaves less than a factor 3 above 4 µm, and the broken power law does not describe its w²(ℓ). Its free exponent is still 0.93, not 0.63.
     - Study 3: short wall, ζ ≈ 0.9.
2. **KPZ term.** Λ is of order one or smaller in every study (≈ 0–1.4): between the elastic (0) and the isotropic (1) values, never ≫ 1.
   - **R stays at ≤ 0.2–0.3:** at scales below l₀ the tilts change the velocity by at most 20–30 %. Even if λ were finite, the nonlinear term would be a perturbation on these scales.
   - **Λ does not grow as 1/v.** Fixed temperature and disorder (Studies 5, 7): v drops 4.4× from 0.42 A to 0.36 A, and Λ goes from ≈ 0.3–0.7 to ≈ 1.2. With a fixed λ (qKPZ) it would grow 4.4×.
   - **Across all studies:** v spans a factor 20 (0.12–2.36 nm/s, three samples) and Λ stays O(1) without a trend that would point to λ/v diverging.
   - The isotropic value Λ = 1 (normal growth, λ = v) is compatible within the spread.
3. **Second sample (P7925, Studies 8 and 9).** Same picture as Pt/Co/Pt.
   - **Study 8:** ζ ≈ 1.05 below l₀ ≈ 14–23 µm, and Λ ≈ 1.
   - **Study 9** is the least clean case:
     - ζ ≈ 0.9 (folds excluded) without a clear break inside the 100 µm field;
     - a linear term in v_n(φ) (b ≈ +0.3), i.e. a preferred direction of motion other than the mean normal;
     - but Λ is still O(1) and R ≤ 0.37.
   - The pixel size of these runs is not known. Exponents, Λ and R do not depend on it.
4. **Third sample, ferrimagnetic GdFeCo (Study 10).** Same picture again.
   - ζ ≈ 1.0–1.2 below l₀ ≈ 18–35 µm, ζ₁ = 1.31 with the crossover fitted; χ² 13 (qEW) against 1266 (qKPZ).
   - Λ ≈ 0.1–0.75 for L = 2–8 µm, R ≤ 0.24. At L = 16 µm the contour gives Λ ≈ +2 and h(x) Λ ≈ −1.4: opposite signs on few segments, not a KPZ signal.
   - It is the fastest wall (2.4 nm/s), and still no Λ ≫ 1.
5. **Study 11 is the exception: ζ ≈ 0.7.** Same sample and drive as Studies 8–9, but the wall advances only 3.7 µm in 6.9 h (mostly in one bulge), about the dynamical correlation length ξ ≈ 2–3 µm. At larger scales the shape is still that of the initial wall, not a driven steady state, and a wall prepared at rest is expected to be in the equilibrium regime (ζ = 2/3); ζ also drifts over the run (0.92 → 0.52). Λ = O(1) and R ≤ 0.17, as everywhere else, so test B gives no sign of qKPZ. Its roughness is not used for the decision.
6. **Conclusion.** Below l₀ the walls behave as qEW interfaces, with at most the geometric KPZ term of an isotropic medium (λ ∝ v). The qKPZ signatures are absent: no ζ ≈ 0.63, and no λ/v that grows as the creep velocity falls.

## Which crossover: ξ_av or l₀? (step11, 2026-10-15)
ζ ≈ 1/2 is expected both above l₀ (effective: the wall is multivalued) and above the creep avalanche scale ξ_av (genuine: the wall is single-valued, thermal regime). The break l× found by step8 is therefore min(ξ_av, l₀), and the exponent cannot tell which. `step11_overhangs.py` measures directly whether the wall is single-valued at each scale: P_ov(ℓ) is the fraction of contour segments of extent ℓ, each in its own frame, that go back by more than 3 px (README, step 11).

| study | l× (step8, local width) [µm] | P_ov(l×) | l₀^geom: P_ov = 1/2 [µm] | P_ov at ½ wall | mean spacing of overhangs [µm] | reading |
|---|---|---|---|---|---|---|
| 5 | 80 | 0.02 | > 97 | 0.01 | 6300 | single-valued at all scales |
| 6 | 13 | 0.04 | 68 (65–70) * | 0.21 | 160 | break ≪ l₀^geom |
| 7 | 15 | 0.02 | > 66 | 0.04 | 500 | break ≪ l₀^geom |
| 2 | 28 | **0.40** | **36 (31–42)** | 0.68 | 33 | **break = overhang crossover** |
| 3 | 10 | 0.21 | 43 (39–46) | 0.39 | 31 | break 4× below l₀^geom |
| 4 | 16 | 0.06 | > 78 | 0.13 | 150 | break ≪ l₀^geom |
| 8 | 23 | 0.01 | > 110 | 0.02 | 410 | break ≪ l₀^geom |
| 9 | 48 | **0.40** | **58 (41–58)** * | 0.41 | 53 | **break ≈ overhang crossover** |
| 10 | 35 | 0.18 | > 51 | 0.14 | 200 | break < l₀^geom |
| 11 | 26 | 0.21 | 46 (44–59) | 0.33 | 72 | no clean break |

\* reached only with windows longer than half the wall, where one large fold can dominate.

- l₀^geom does not depend on the tolerance (2, 3 or 5 px change it by ≤ 15 %).
- **Only in Studies 2 and 9** is the wall about half multivalued at the break (P_ov ≈ 0.4). These are, with Study 3, the runs where overhangs are common (one every 30–50 µm along the wall). There the break is the overhang crossover of Kolton, Ferrero & Rosso.
- **In Studies 4–8 and 10** the wall is still single-valued at the break (P_ov ≤ 0.2, mostly ≤ 0.06), and overhangs are rare (one every 150–6000 µm). The break there is not l₀; l₀ lies above it, often beyond the field of view.
- A break in a single-valued wall is what the avalanche scale ξ_av would produce (thermal regime, ζ = 1/2 above it). It is a candidate, not a proof. Other causes of a break in a single-valued wall: the macroscopic shape of the wall in the field of view (curvature set by the nucleation and the aperture), and the limited range of ℓ near the wall length.
- Studies 5–7 (same sample and T): l× = 80, 13, 15 µm at 20.7, 21.9, 24.2 Oe, all on single-valued walls. If l× is ξ_av, it falls with the field, the expected direction, but much more steeply than L_opt ∝ H^−3/4. The test that would settle it is several fields at fixed T in the same region of the film.
- **The qEW/qKPZ tests are not affected.** The window [4 µm, l×] is depinning-like whichever crossover closes it, and in most studies the wall is single-valued in the whole window.

## Bubbles and the scale of the fit (Studies 12–14, 2026-10-16)
Expanding bubbles are often used for this question; a recent one (Durin et al., PRB 110, L020405 (2024): CoFeB, 0.13–0.16 mT, 0.4 µm/px) reads ζ ≈ 0.63 from the polar S(q) and from cluster shapes as qKPZ. Studies 12–14 (group movies, `arrival_map.py` + `contour_analysis.py`) give:

| study | sample, drive, resolution | scales of the fit | ζ local width, raw → pixel-bias corrected | P_ov | Λ local (step10 B) | Λ(φ) macroscopic |
|---|---|---|---|---|---|---|
| 12 | Pt/(CoNi)×3 bubble, 10 ms pulses, 0.467 µm/px | 3.7–37 µm | 0.75 → **≈ 0.64** (20–80 px: 0.69 → 0.67) | 0 | 0.43–1.36 | 0.86–1.19 |
| 13 | PtCoNiAl bubble, 150 Oe, 0.467 µm/px, 9 frames | 3.7–25 µm | 0.68 → ≈ 0.5 | 0 | 1.0–1.07 | 0.83–1.15 |
| 14 | PtCoNiAl, two creep walls, 10.4 Oe, 0.117 µm/px | 0.9–9 µm | 0.91 → ≈ 0.7–0.85 | ≤ 0.15 | 0.79–1.08 | – |

- The bias correction comes from synthetic bubbles of known ζ and the same amplitude analysed with the same code (`synthetic_bubble.py`, table in the Study 12 report): walls only ~1 px rough make the local-width exponent too large for ζ < 1 and too small for 1.25, and the bias grows with the amplitude.
- **A bubble in a Pt/(CoNi) multilayer gives ζ ≈ 0.64 with isotropic growth and Λ = O(1)**, i.e. the bubble exponent alone is not evidence for qKPZ. With isotropic growth the KPZ coefficient is the geometric one, λ ≈ v (Λ ≈ 1), which vanishes at depinning; qKPZ needs an orientation-dependent depinning threshold, and in creep (velocity exponentially sensitive to the barrier) that would make the bubble facet.
- **The exponent depends on the scales of the fit:** ≈ 0.5–0.64 for the bubbles over tens of µm at 0.47 µm/px, 0.7–0.9 for the walls of Study 14 over 1–9 µm, 0.9–1.25 for Studies 1–11 over 4 µm–l₀. That is what a crossover from a depinning-like exponent at small scales to ζ = 1/2 above ξ_av or l₀ produces when fits at different resolutions straddle it. The samples and drives differ, so this is a consistent picture, not a test; the test is a bubble at 0.12 µm/px, or a bubble series in field.

## Caveats
- **v_n(φ) is noisy for the slowest walls.** In Studies 3 and 6 the wall moves ~4–8 µm in 10–11 h, so v_n(φ) rests on few switching events. Some Λ values there are negative (velocity falling faster than cos φ with the tilt), with wide errors.
- **Pinned segments.** Segments pinned at defects are tilted and slow. With the defect mask Λ changes by ~0.1–0.6 (Study 2: 0.6–1.2 masked vs 0.3–0.8 unmasked) but stays O(1).
- **The tests are on the creep motion.** Within an avalanche the dynamics is depinning-like, but the 20 s frames do not resolve it. B measures the coarse-grained velocity over 2–15 px of advance.
- **Not used for the decision:**
  - test C, the avalanche shape (mean advance S/ℓ_y ∝ ℓ_y^0.3–0.7): the avalanches are slivers at the resolution;
  - test D, the skewness of the local heights (≈ −0.1 to −0.4 at 20 µm, lagging valleys): pinning tongues also produce it.
- **The step8 crossover differs between runs of the same sample:** ≈ 70–80 µm at 0.36 A against 10–30 µm in the others. step11 shows that in Studies 5–7 the wall is single-valued at all these breaks, so they are not l₀ (see above); a field-dependent ξ_av, or a different region of the film (Study 5 was recorded two weeks after Studies 6 and 7), are both possible.
