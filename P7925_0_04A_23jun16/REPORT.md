# Study 11: `P7925_0_04A_23jun16` (sample P7925, 0.04 A, Δt = 8.25 s, 6.9 h)

First report, 2026-10-15. Pipeline and definitions: see the [README](../README.md). Same sample and drive as [Study 8](../P7925_0_04A_24jun16/REPORT.md) (one day later) and [Study 9](../P7925_0_04A_26jun16/REPORT.md). Independent of the other studies.

## Data

- **Source:** `/media/…/MOBILEDATAMOKE/DataDWMoviesLargas/23jun16/`, copied to `~/MOKEdata_local/23jun16` (not in the repository).
  - 3000 TIFF frames, 1392×1040 px, one wall.
  - `info.txt`: I = 0.04 A, exposure 2.6 s.
- **Δt = 8.25 s:** from the file modification times (uniform, no pauses; `--fps 0.12121`).
- **Scale 0.1175 µm/px: assumed**, as in Studies 8 and 9.
- **Geometry:** as Study 8, horizontal wall, bright domain growing upward: `--rotate 3 --bright-domain`.

```bash
python step1_extract_wall.py ~/MOKEdata_local/23jun16 --outdir P7925_0_04A_23jun16 --fps 0.12121 --um-per-px 0.1175 \
       --rotate 3 --bright-domain --offset auto --flatfield 2 --init-wall --swept-thr 0.018 --rows 130 980
python step2_persistence.py P7925_0_04A_23jun16 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py P7925_0_04A_23jun16 --every 10 --rmin 25
python step8_crossover.py P7925_0_04A_23jun16
python step10_universality.py P7925_0_04A_23jun16
python step11_overhangs.py P7925_0_04A_23jun16
```

## Wall extraction and checks

**Segmentation check.** `fig_qc_overlay.png` and `qc/overlay_*.png` (1 frame in 50) were approved by the user before the analysis.

Two options were needed on top of those of Study 8:
- **`--swept-thr 0.018`.** The wall moves so little that the swept band is a small minority of the field. The switching contrast has a noise peak at 0 and the swept band at ≈ 0.03, and Otsu put the threshold at 0.0023, inside the noise peak: 38 % of the aperture was marked as swept and the wall escaped to the aperture edge after frame ~1500. 0.018 is the minimum between the two peaks.
- **`--rows 130 980`.** Below row ~985 (bottom corner of the aperture) the contrast is too weak and the detected wall leaves the visible one.

- **Rows:** 130–979 (850 rows = 100 µm).
- **Drift:** dy = 0…+3.2 px, dx = 0…+0.8 px.
- **Motion:** all 3000 frames used. Total advance 32 px = 3.7 µm, **mean velocity 0.15 nm/s**, half that of Study 8; no retreats > 1 px. The advance is very uneven along the wall (`fig_kymograph.png`): the middle (rows ~250–550) bulges forward by up to ~100 px, the bottom part barely moves.
- **Noise:** σ_Δh = 0.52 px.
- **Tilt:** θ0 = −7.4° in step 7; overhangs on the base plane 1.6 % of (row, frame).

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted. τ in frames of 8.25 s.

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] at τ |
|---|---|---|---|---|
| 0.5 | 370 (51 min) | 9.7 | 1.74 ± 0.08 | 2.45 at 2.0 h |
| 1 | 370 (51 min) | 10.3 | 1.93 ± 0.08 | 2.76 at 2.0 h |
| **1.5** | **355 (49 min)** | 14.4 | **2.28 ± 0.09** | 3.28 at 2.0 h |
| 2 | 370 (51 min) | 14.0 | 2.20 ± 0.08 | 3.25 at 2.0 h |
| 3 | 413 (57 min) | 13.1 | 2.08 ± 0.08 | 3.22 at 2.0 h |
| 5 | 475 (65 min) | 11.5 | 2.10 ± 0.08 | 3.17 at 2.2 h |

ξ(τ\*) ≈ 2.1–2.3 µm and τ\* ≈ 50–60 min, as in Study 8 (2.3–2.5 µm, 35–50 min).

## Roughness and l₀

`step7` (`fig_roughness.png`, `fig_local_width.png`): 852 rows on the base plane (100 µm), 300 frames, fit range 2.9–25 µm.

| estimator | ζ |
|---|---|
| global S(q), u_area | 0.68 ± 0.05 (thirds: 0.92, 0.66, 0.52) |
| local width, local rotation (contour) | 0.69 ± 0.02 (folds excluded 0.78, tilts > 45° excluded 0.71) |
| local width on u_area | 0.73 |
| local S(q, l), l = 14–55 µm | 0.65–0.84 |

- **ζ ≈ 0.7 from 1 to ~30 µm**, a single power law. Above ~30 µm w²(ℓ) steepens: the macroscopic bulge of the wall.
- **No break to ζ ≈ 0.5.** step8: l₀ = 26 µm from the local width (variants 13–90 µm), 4.5 µm from the global S(q), 5 µm from the local slope. None is a clean crossover.
- **Height distribution:** skewness −1.2 on the base plane (the bulge), ≈ 0 in 5–20 µm windows.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A. Roughness** (ζ over [4 µm, l₀ = 26 µm]):

| estimator | ζ |
|---|---|
| contour | 0.69 (0.67–0.71) |
| folds excluded | 0.77 (0.76–0.79) |
| u_area | 0.71 (0.70–0.73) |

- Broken power law (folds excluded): ζ₁ free = 0.89, with l₀ pushed to the end of the field (no break). **χ² qEW / qKPZ / equilibrium = 1700 / 3714 / 2827**, all poor fits (free: 384).

**B. KPZ term:**
- Λ from the contour: 0.3–1.2 for L = 2–8 µm (median 0.65); at L = 16 µm the intervals span −18 to +12 (unusable).
- Λ from h(x): −0.4 to 0.6 (median ≈ 0.1).
- **R ≤ 0.17** where Λ is measured.

**Overhangs (step11):** P_ov = 0.05 at 5 µm, 0.17 at 20 µm; l₀^geom = 46 µm (44–59). At the step8 break (26 µm) P_ov = 0.21.

## Reading

This is the only study with ζ ≈ 0.7 below l₀, between the creep equilibrium value 2/3 and qKPZ 0.63, instead of 1.0–1.25. The KPZ test does not point to qKPZ: Λ = O(1) and R ≤ 0.17, as in all the other studies.

The likely cause is that **the wall hardly moves**. It advances 3.7 µm on average in 6.9 h, mostly in one bulge, while the correlation length of the dynamics is ξ ≈ 2–3 µm. At scales of more than a few µm the measured shape is therefore mostly that of the initial wall, as it was prepared before the run, not the steady state of the driven wall. A wall that relaxed at rest (or in very slow creep) is expected to be in the equilibrium regime, ζ = 2/3, at those scales. The drift of ζ over the run (thirds: 0.92, 0.66, 0.52) also says that the shape is not stationary.

**This run should not be used for the qEW/qKPZ decision on its roughness**, only for the KPZ term at small scales (which agrees with the other studies).

## Caveats

- **Scale assumed** (0.1175 µm/px); exponents, Λ and R do not depend on it.
- **Manual threshold** (`--swept-thr 0.018`) and **rows restricted** to 130–979.
- **Not run yet:** steps 3–6 and 9.
