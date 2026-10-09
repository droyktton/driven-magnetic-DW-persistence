# Study 3: `PtCoPt_0_48A_1sep16` (Pt/Co/Pt, 0.48 A, Δt = 22.05 s, 11 h)

First report, 2026-10-09. Pipeline and definitions: see the [README](../README.md). Same sample, microscope and protocol as [Study 4](../PtCoPt_0_50A_2sep16/REPORT.md) (0.50 A, the next day), which explains the acquisition and the step 1 options in more detail. Independent of Studies 1 and 2.

## Data

- Folder `/media/…/MOBILEDATAMOKE/MOKE 2016/1sep16/` (not in the repository): 1800 TIFF frames, 1392×1040 px, one wall driven at 0.48 A.
- **Δt = 22.05 s** from the file modification times; 11.0 h. **Scale 0.1175 µm/px** (20x, confirmed by the user). Camera offset ≈ 33 000 counts.
- **The third run of the series, 30ago16 (0.42 A), is not analysed.** Its frame interval is 22.25 s with a ~50 s pause every 100 frames (mean 22.75 s). The wall barely moves in 11 h, and the absolute contrast between the domains changes and fades during the run. So the bright/dark reference method finds the domain pattern instead of a swept band. It would need a frame-by-frame detector of the intensity step.

```bash
python step1_extract_wall.py "/media/…/MOKE 2016/1sep16" --outdir PtCoPt_0_48A_1sep16 --fps 0.045351 \
       --um-per-px 0.1175 --offset auto --flatfield 2 --crop-left 400 --init-wall --rows 260 840
python step2_persistence.py PtCoPt_0_48A_1sep16 --row-mean --eps 0.5 1 1.5 2 3 5
python step3_eps_sweep.py PtCoPt_0_48A_1sep16 --row-mean
python step4_avalanches.py PtCoPt_0_48A_1sep16 --tau-m 1 3 10 30
python step5_pinning_map.py PtCoPt_0_48A_1sep16 --eps 1.5 --tau 209
python step7_roughness.py PtCoPt_0_48A_1sep16 --every 5 --rmin 25
python step8_crossover.py PtCoPt_0_50A_2sep16 PtCoPt_0_48A_1sep16
python step9_wavelet.py PtCoPt_0_48A_1sep16
```

## Wall extraction and checks

- **Rows 260–839 only** (580 rows, 68 µm). In the unrestricted run, above y ≈ 250 a spurious contrast patch appears after frame ~1000 and the line drifts right of the boundary. Below y ≈ 850 the wall jumps along defects late in the run. Checked with the user.
- **Segmentation check:** `fig_qc_overlay.png`, `qc/overlay.mp4`, `qc/overlay_NNNN.png`. Approved by the user. The last frames show small wiggles near y ≈ 300–380.
- **Drift:** dy = −15.0…+0.5 px, dx = −0.6…+2.5 px, with ±5 px jumps over a few tens of frames after frame ~1000. The registration follows them (scatter 0.20 / 0.05 px), but this is the least reliable part of the run.
- **Motion:** all 1800 frames used. Total advance 67 px = 7.8 µm, **mean velocity 0.04 px/frame = 0.21 nm/s**, half that of Study 4 and 7× slower than Study 2.
- **Noise:** σ_Δh = 0.57 px, so 3σ_Δh = 1.7 px. In this run ε ≥ 1.5 px is above the noise.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted. τ in frames (22.05 s). ξ_max is taken over τ < 600. Beyond that, Π < 0.2 and ξ(τ) rises again, which is noise at the edge of the τ range (`fig_xi_chi4_sub_rm_eps1.5.png`).

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | ξ_max [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.85 | 62 | 26 | 11.2 | 11.6 ± 0.4 (1.37) | 14.5 (1.70) at 99 |
| 1 (0.12) | 0.92 | 85 | 65 | 12.4 | 14.2 ± 0.5 (1.67) | 14.7 (1.73) at 93 |
| **1.5 (0.18)** | 0.99 | 219 | **209 (77 min)** | 29.7 | **21.4 ± 0.8 (2.52)** | 21.5 (2.53) at 204 |
| 2 (0.24) | 0.99 | 227 | 213 (78 min) | 29.9 | 21.4 ± 0.8 (2.52) | 21.5 (2.53) at 204 |
| 3 (0.35) | 1.00 | 257 | 221 (81 min) | 30.3 | 21.6 ± 0.9 (2.54) | 21.8 (2.56) at 207 |
| 5 (0.59) | 1.00 | 330 | 292 (107 min) | 31.4 | 21.3 ± 0.8 (2.50) | 21.9 (2.57) at 213 |

- **Above the noise, τ\* ≈ 77 min** (2× Study 4, consistent with half the velocity) and **ξ(τ\*) ≈ 2.5 µm**, independent of ε.
- **Here ξ(τ) peaks at τ\***, so ξ(τ\*) ≈ ξ_max. ξ(τ\*) is the same as in Study 4 (2.6 µm). ξ_max is smaller than Study 4's (3.0 µm).
- Between ε = 1 and 1.5 px everything jumps (τ\* 65 → 209 frames, χ4 12 → 30): below 3σ_Δh the persistence is set by the noise.

## Avalanches

| τ_m | avalanches ≥ 20 px | share of swept area | median S [µm²] | max S [µm²] | τ (with cutoff) | S_cut [µm²] |
|---|---|---|---|---|---|---|
| 22 s | 373 | 42 % | 0.48 | 3.4 | 0.90 ± 0.31 | 0.55 |
| 66 s | 421 | 65 % | 0.59 | 5.3 | 0.97 ± 0.19 | 1.1 |
| 220 s | 340 | 82 % | 0.83 | 14.9 | 0.87 ± 0.14 | 2.2 |
| 660 s | 231 | 92 % | 0.99 | 26.8 | 1.16 ± 0.11 | 7.6 |

At 22 s only 42 % of the swept area is in resolved avalanches; the majority advances in sub-resolution steps. The distributions are again cutoff-dominated.

## Pinning map

`step5_pinning_map.py --eps 1.5 --tau 209`:

- Velocity in laguna pixels is 0.53 × that in avalanche pixels.
- Per-row persistence m_i = 0.51 ± 0.21, strongly anticorrelated with the row velocity (r = −0.70). Its correlation length along the wall is 56 px = 6.6 µm.
- Rows with more lagunas are *less* persistent (r = −0.28), as in Study 2. In Study 4 the sign is opposite (+0.23).

## Roughness

Base plane θ0 = 18.7°, 593 rows (70 µm), 360 frames; overhangs 4.0 % of row-frames. Fit range 2.9–17 µm.

| estimator | ζ |
|---|---|
| global S(q), u_area | 0.90 ± 0.11 (thirds: 0.97, 0.82, 0.91) |
| global S(q), u_front / u_back | 0.91 ± 0.10 / 0.78 ± 0.11 |
| local width, local rotation (contour) | 0.69 ± 0.03 (folds excluded 0.91; u_area 0.88) |
| local S(q, ℓ), ℓ = 13–40 µm | 0.75–1.0 |

- **ζ ≈ 0.8–0.9 over 3–17 µm**, between Study 2 (0.75) and Study 4 (1.0). Here the thirds of the run agree better.
- **l₀** (`step8`): 10 µm from the local width (variants 2–110 µm, 4/30 not bracketed), 7.3 µm from the global S(q), 40 µm from the local slope; l_fold = 17 µm, l_tilt = 7 µm. Poorly constrained: the wall is only 70 µm long.
- **Wavelets:** no usable scale (kept fraction ≤ 0.77), no ζ reported.
- **Height distribution:** width 1.37 µm, skewness −0.14 ± 0.10, excess kurtosis −0.4 ± 0.3.

## Comparison of the field series (Studies 3 and 4)

| | Study 3 (0.48 A) | Study 4 (0.50 A) | Study 2 (0.8 A, 2015, other conditions) |
|---|---|---|---|
| v [nm/s] | 0.21 | 0.43 | 1.65 |
| τ\* (ε ≈ 3σ_Δh) | 77 min | 41 min | 2.3 min |
| ξ(τ\*) [µm] | 2.5 | 2.6 | 1.6 |
| ξ_max [µm] | 2.5 | 3.0 | 2.6 |
| ζ, global S(q), 3–20 µm | 0.90 ± 0.11 | 1.03 ± 0.11 | 0.82 ± 0.09 |
| l₀, local width [µm] | ~10 (poor) | 15 | 16–28 |

A 4 % change in drive current doubles the velocity and halves τ\*, as expected for creep. ξ(τ\*) does not change. The two runs image the same field of view (the same defects appear in both), but each day starts from a new wall: the wall of 2sep16 starts to the left of where that of 1sep16 ended. With two points no field dependence of ξ can be claimed.
