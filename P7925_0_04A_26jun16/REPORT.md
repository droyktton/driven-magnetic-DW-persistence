# Study 9: `P7925_0_04A_26jun16` (sample P7925, 0.04 A, Δt ≈ 7.8 s, 13 h)

First report, 2026-10-14. Pipeline and definitions: see the [README](../README.md). Same sample, drive and microscope as [Study 8](../P7925_0_04A_24jun16/REPORT.md), two days later; the two are compared at the end.

## Data

- **Source:** folder `/media/…/MOBILEDATAMOKE/DataDWMoviesLargas/26jun16/`, copied to `~/MOKEdata_local/26jun16` (not in the repository). 6000 TIFF frames, 1392×1040 px, I = 0.04 A (`info.txt`).
- **Δt:** median 7.25 s, mean 7.79 s, with a pause every 150 frames. The mean is used (`--fps 0.12844`).
- **Scale 0.1175 µm/px: assumed**, as in Study 8.
- **Geometry:** the wall is horizontal, and here the dark domain grows downward. Used `--rotate 1` without `--bright-domain`.

```bash
python step1_extract_wall.py ~/MOKEdata_local/26jun16 --outdir P7925_0_04A_26jun16 --fps 0.12844 --um-per-px 0.1175 \
       --rotate 1 --offset auto --flatfield 2 --init-wall --rows 0 1150
python step2_persistence.py P7925_0_04A_26jun16 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py P7925_0_04A_26jun16 --every 20 --rmin 25
python step8_crossover.py P7925_0_04A_24jun16 P7925_0_04A_26jun16
python step10_universality.py P7925_0_04A_24jun16 P7925_0_04A_26jun16
```

## Wall extraction and checks

- **Initial domain.** The hysteresis threshold of the swept area leaked into the initial domain (left of the `--init-wall` path). This is the reason for the new step 1 rule: the initial domain cannot switch.
- **`--rows 0 1150`.** The rows below are near the edge of the aperture. The output keeps image rows 306–1149 (844 rows = 99 µm).
- **Segmentation check.** `fig_qc_overlay.png` and `qc/overlay_small.mp4` were approved by the user before the analysis.
- **Drift:** dy = −4.3…+0.7 px, dx = −0.2…+1.0 px. The drift jumps by 3–4 px when the focus is readjusted (`info.txt` lists a focus-voltage table). Scatter around the smoothed curve 0.10 / 0.02 px.
- **Motion:** all 6000 frames used. Total advance 76 px = 8.9 µm, **mean velocity 0.013 px/frame ≈ 0.19 nm/s** (0.63× Study 8), no retreats > 1 px.
- **Noise:** σ_Δh = 0.57 px.
- **Tilt:** θ ≈ −15° to −17°.
- **Overhangs:** 2.9 % of (row, frame), 10× Study 8.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted. τ in frames (≈ 7.8 s).

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] at τ |
|---|---|---|---|---|
| 0.5 | 292 (38 min) | 11.1 | 1.80 ± 0.07 | 4.70 at 4.8 h |
| 1 | 418 (54 min) | 11.7 | 2.21 ± 0.06 | 4.75 at 4.8 h |
| **1.5** | **243 (32 min)** | 13.0 | **1.96 ± 0.03** | 4.97 at 4.8 h |
| 2 | 242 (31 min) | 12.1 | 1.87 ± 0.02 | 4.88 at 4.8 h |
| 3 | 338 (44 min) | 11.1 | 2.01 ± 0.03 | 4.42 at 5.1 h |
| 5 | 518 (67 min) | 10.2 | 2.25 ± 0.04 | 3.60 at 5.4 h |

- **ξ(τ\*) ≈ 1.9–2.3 µm**, as in Study 8.
- **ξ_max ≈ 4.4–5.0 µm**, larger than in any other study. It is reached at τ ≈ 4.8 h, where few sites are still persistent.
- τ\* does not grow monotonically with ε (243 frames at 1.5 px against 418 at 1 px): χ4(τ) has a broad maximum.

## Roughness and l₀

`step7`: 855 rows on the base plane (100 µm), 300 frames, fit range 2.9–25 µm.

| estimator | ζ |
|---|---|
| global S(q), u_area | 0.74 ± 0.03 (thirds: 0.76, 0.76, 0.67) |
| local width, local rotation (contour) | 0.83 ± 0.03 (folds excluded 0.96, tilts > 45° excluded 0.92) |
| u_area | 0.89 |
| B(r) | 0.76 |

- **Folds lower ζ.** Here the folds lower the contour exponent much more than in Study 8. Excluding them (9 % of the segments) raises ζ from 0.83 to 0.96.
- **Crossover l₀** (`step8`):

  | estimator | l₀ [µm] | bootstrap | variants |
  |---|---|---|---|
  | local width | 48 | 44–52 | 26–127, 6/30 not bracketed |
  | global S(q) | 29 | — | 27–128 |
  | local slope | 0.8 | — | — |

  - The local slope gives no useful crossover.
  - Proxies: l_fold = 55 µm, l_tilt = 38 µm.
  - **l₀ ≈ 30–50 µm, poorly bounded.** w²(ℓ) keeps growing with ζ ≈ 1 up to ~100 µm, close to the size of the field.
- **Height distribution** (20 µm windows): skewness −0.39 ± 0.10, excess kurtosis +0.1 ± 0.3.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A. Roughness below l₀** (ζ over [4 µm, l₀ = 48 µm]):

| estimator | ζ |
|---|---|
| contour | 0.81 (0.78–0.85) |
| folds excluded | 0.90 (0.89–0.91); 0.94 up to the lower bound of l₀ |
| u_area | 0.88 |

- **Broken power law:** ζ₁ free ≈ 1.0 (pca 1.03, folds excluded 1.01, u_area 0.93), but with l₀ at the edge of the allowed range (400 µm): no break inside the field.
- **χ² (folds excluded):**

  | model | χ² |
  |---|---|
  | qEW | 1161 |
  | qKPZ | 3948 |
  | equilibrium | 3312 |
  | free | 459 |

  qEW is still preferred over qKPZ, but it fits worse than the free exponent.

**B. KPZ term:**
- **Λ:**

  | segments | Λ |
  |---|---|
  | contour, defects masked | 0.7–2.0 (median 1.2) |
  | contour, unmasked | 0.4–1.1 |
  | h(x) | 0.5–0.75 |

- **R ≤ 0.37.**
- v_n(φ) has a clear linear term b ≈ +0.3. Segments tilted one way move faster than the others, which suggests that the wall moves partly along a direction other than its mean normal (tilted drive, or the aperture edge). The quadratic coefficient Λ is O(1).

**Conclusion: closer to qEW than to qKPZ.**
- Λ is O(1) and does not grow, and ζ ≈ 0.9–1.0 is far from 0.63.
- But ζ is lower than 1.25 and l₀ is not well defined: this run is less conclusive than Study 8.

## Comparison with Study 8

| | Study 8 (24jun16) | Study 9 (26jun16) |
|---|---|---|
| v [nm/s] | 0.30 | 0.19 |
| ξ(τ\*) / ξ_max [µm] | 2.3 / 2.9 | 2.0 / 5.0 |
| overhangs | 0.3 % | 2.9 % |
| ζ below l₀, folds excluded | 1.06 | 0.90 |
| l₀ [µm] | 14–23 | 30–50 (poorly bounded) |
| Λ contour / h(x) | 1.4 / 0.8 | 1.2 / 0.6 |
| max R | 0.25 | 0.37 |

Same sample and nominal drive, but the wall of Study 9 is slower, more folded and rougher at large scale. Each run starts from a new wall, possibly in a different region of the sample.

## Caveats

- **Assumed scale** and non-uniform Δt, as in Study 8.
- **Focus jumps** of 3–4 px in the drift. They are registered, but frames near a jump may be slightly defocused.
- **Not run yet:** steps 3–6 and 9.
