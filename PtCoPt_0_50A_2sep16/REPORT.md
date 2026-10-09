# Study 4: `PtCoPt_0_50A_2sep16` (Pt/Co/Pt, 0.50 A, Δt = 22.05 s, 11 h)

First report, 2026-10-09. Pipeline and definitions: see the [README](../README.md). Same sample, microscope and protocol as [Study 3](../PtCoPt_0_48A_1sep16/REPORT.md) (0.48 A, the day before); the two are compared at the end. Independent of Studies 1 and 2.

## Data

- Folder `/media/…/MOBILEDATAMOKE/MOKE 2016/2sep16/` (not in the repository): 1799 TIFF frames, 1392×1040 px, of one wall driven at 0.50 A. Part of a field series with 30ago16 (0.42 A, not analysed, see Study 3) and 1sep16 (0.48 A, Study 3).
- **Δt = 22.05 s** from the file modification times (no Micro-Manager metadata). Total duration 11.0 h.
- **Scale 0.1175 µm/px** (20x, the same microscope as Study 2), confirmed by the user; the files carry no pixel size.
- **Camera offset ≈ 33 000 counts**, against a signal of only ~7000 above it. It must be subtracted before any ratio.
- An octagonal illuminated aperture fixed on the camera; its illumination profile changes during the run. Many static round defects.

```bash
python step1_extract_wall.py "/media/…/MOKE 2016/2sep16" --outdir PtCoPt_0_50A_2sep16 --fps 0.045351 \
       --um-per-px 0.1175 --offset auto --flatfield 2 --crop-left 400 --init-wall --rows 0 815
python step2_persistence.py PtCoPt_0_50A_2sep16 --row-mean --eps 0.5 1 1.5 2 3 5
python step3_eps_sweep.py PtCoPt_0_50A_2sep16 --row-mean
python step4_avalanches.py PtCoPt_0_50A_2sep16 --tau-m 1 3 10 30
python step5_pinning_map.py PtCoPt_0_50A_2sep16 --eps 1.5 --tau 112
python step7_roughness.py PtCoPt_0_50A_2sep16 --every 5 --rmin 25
python step8_crossover.py PtCoPt_0_50A_2sep16 PtCoPt_0_48A_1sep16
python step9_wavelet.py PtCoPt_0_50A_2sep16
```

## Wall extraction and checks

The new step 1 options are described in the README (*Options for aperture images with a camera offset*). Why each one was needed here:

- `--offset auto`: without it the contrast map and the frame selection were wrong (only frames 1089–1799 were kept).
- Aperture mask (automatic with `--offset`): otherwise the Otsu threshold picks the moving edge of the aperture.
- `--flatfield 2`: the illumination drifts non-uniformly; the vignetted left side showed fake switching contrast.
- `--crop-left 400`: the polynomial does not follow the aperture edge on the far left.
- `--init-wall`: parts of the wall never move in 11 h; without a directly located initial wall those rows had no wall at all.
- `--rows 0 815`: below y ≈ 820 the wall runs almost horizontally along a line of defects near the aperture edge late in the run (`fig_qc_overlay.png` of the unrestricted run, checked with the user). The output keeps rows 75–814 (740 rows, 87 µm).

**Segmentation check.** `fig_qc_overlay.png`, `qc/overlay.mp4` (all frames) and `qc/overlay_NNNN.png`: the line follows the dark/bright boundary in all frames, including around the defects. Approved by the user before the analysis.

- **Drift:** dy = −10.6…+0.3 px, dx = −0.5…+2.6 px, with fast jitter after frame ~1250; scatter around the smoothed curve 0.21 / 0.11 px (larger than in Study 2, 0.03–0.04 px).
- **Motion:** all 1799 frames used. Total advance 143 px = 16.8 µm, **mean velocity 0.08 px/frame = 0.43 nm/s**, no retreats > 1 px. The advance is not uniform: slower between frames ~600 and 1200, faster afterwards (`fig_h_mean.png`).
- **Noise:** σ_Δh = 0.52 px (backward one-frame steps), so ε = 3σ_Δh ≈ 1.5 px is the reference threshold. The noise is 1.6× that of Study 2 while the velocity is 3.5× smaller: in one frame the wall moves much less than the noise.
- **Tilt:** θ ≈ 17–19° from vertical on the base plane (θ0 = 17.2° in step 7). As in Study 2, h is measured along x without the cos θ correction (≈ 5 % here).

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted (`*_sub_rm_*`). τ in frames (1 frame = 22.05 s).

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | ξ_max [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.81 | 20 | 10 (3.7 min) | 17.2 | 12.9 ± 0.4 (1.51) | 20.3 (2.38) at 175 |
| 1 (0.12) | 0.93 | 37 | 29 (11 min) | 18.8 | 17.2 ± 0.5 (2.02) | 20.1 (2.36) at 88 |
| **1.5 (0.18)** | 0.98 | 65 | **112 (41 min)** | 21.9 | **22.1 ± 1.2 (2.59)** | 25.9 (3.04) at 436 |
| 2 (0.24) | 0.99 | 70 | 112 (41 min) | 22.6 | 22.0 ± 1.2 (2.58) | 25.7 (3.02) at 441 |
| 3 (0.35) | 0.99 | 84 | 122 (45 min) | 24.3 | 21.7 ± 1.3 (2.54) | 25.2 (2.96) at 366 |
| 5 (0.59) | 1.00 | 111 | 143 (53 min) | 27.1 | 21.5 ± 1.4 (2.53) | 24.9 (2.92) at 381 |

Figures: `fig_eps_sweep_sub_rm.png`, `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`.

1. **Above the noise (ε ≥ 3σ_Δh) the results do not depend on ε.** ξ(τ\*) = 2.5–2.6 µm and ξ_max ≈ 3.0 µm for ε = 1.5–5 px. Below 3σ_Δh, τ\* and ξ(τ\*) are set by the noise: Π(τ=1) < 0.95 and τ\* jumps from 29 to 112 frames between ε = 1 and 1.5 px.
2. **τ\* ≈ 41 min**, 18× longer than in Study 2 (2.3 min): the wall is much slower. Here τ\* comes after Π = ½ (65–111 frames), not at it as in Study 2.
3. **ξ(τ\*) ≈ 2.6 µm and ξ_max ≈ 3.0 µm (at τ ≈ 6–7 h)**, against 1.6 and 2.6 µm in Study 2. ξ_max is reached late, when Π is small (few persistent sites), so it has less statistics.

## Avalanches

`step4_avalanches.py` (`fig_avalanches.png`, `fig_avalanche_map.png`), patches ≥ 20 px, defects excluded.

| τ_m | avalanches ≥ 20 px | share of swept area | median S [µm²] | max S [µm²] | τ (with cutoff) | S_cut [µm²] |
|---|---|---|---|---|---|---|
| 22 s | 1278 | 57 % | 0.48 | 4.8 | 1.27 ± 0.13 | 0.89 |
| 66 s | 1200 | 78 % | 0.64 | 7.8 | 1.03 ± 0.09 | 1.6 |
| 220 s | 788 | 92 % | 1.05 | 13.6 | 0.82 ± 0.08 | 3.1 |
| 660 s | 381 | 97 % | 2.18 | 29.7 | 0.68 ± 0.08 | 7.3 |

As in Study 2, P(S) is dominated by its cutoff (S_cut only 2–4× the smallest resolved area at short windows) and S_cut grows with the window: the fitted τ is not an avalanche exponent. At 22 s, 43 % of the swept area advances in sub-resolution steps (lagunas), more than in Study 2 (34 %).

## Pinning map

`step5_pinning_map.py --eps 1.5 --tau 112` (`fig_pinning_map.png`).

- Local velocity in laguna pixels is 0.66 × that in avalanche pixels (Study 2: 0.72 at the same 4 px smoothing): the lagunas again sit on slow regions.
- Per-row persistence m_i = 0.37 ± 0.14, strongly anticorrelated with the row velocity (r = −0.62; Study 2: −0.41). Its correlation length along the wall is 43 px = 5.1 µm. Laguna fraction per row: correlation length 1.7 µm.

## Roughness

`step7_roughness.py` (`fig_roughness.png`, `fig_local_width.png`, `fig_local_sq.png`, `fig_height_distribution.png`). Base plane θ0 = 17.2°, 741 rows (87 µm), 360 frames. Overhangs are rare here (0.8 % of row-frames, Study 2: 4.7 %).

| estimator (2.9–22 µm) | ζ |
|---|---|
| global S(q), u_area | 1.03 ± 0.11 (thirds of the run: 1.11, 1.30, 0.64) |
| global S(q), u_front / u_back | 1.08 ± 0.09 / 1.08 ± 0.09 |
| local width, local rotation (contour) | 1.03 ± 0.04 (folds excluded 1.11) |
| local S(q, ℓ), ℓ = 13–40 µm | 0.9–1.3 |
| B(r) | 0.40 (biased low for large ζ, see the README self-test) |

- **ζ ≈ 1.0 over 3–22 µm**, larger than Study 2's 0.75 over the same range, between the equilibrium (2/3) and depinning (1.25) values. It fluctuates between thirds of the run (0.64–1.30), as in Study 2.
- **Crossover l₀** (`step8`, `fig_crossover.png`): 15.5 µm from the local width (bootstrap 14–16, variants 11–25), 5.2 µm from the global S(q), 18 µm from the local slope; proxies l_fold = 10 µm, l_tilt = 16 µm. So l₀ ≈ 5–20 µm, at the low end of Study 2 (10–30 µm).
- **Wavelets** (`step9`, `fig_wavelet.png`): no scale passes the criteria (kept fraction < 0.92 at all scales), so no wavelet ζ is reported, as in Study 2. The raw local slopes of F(a) are 1.0–1.5.
- **Height distribution:** width 2.05 µm, skewness −0.35 ± 0.16, excess kurtosis −0.2 ± 0.4 (close to Gaussian).

## Caveats

- Not run yet: `step6_pinning_tests.py` (masks around the defects) and the integer-h cross-check.
- The bottom of the wall (y > 815) is excluded. The top rows are fine.
- Frame-to-frame motion is below the noise: everything relies on lags of tens to hundreds of frames. ξ_max is measured at Π ≈ 0.1–0.2.
