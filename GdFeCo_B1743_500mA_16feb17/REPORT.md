# Study 17: GdFeCo B1743, 500 mA, 16 Feb 2017

Steps 1, 2, 7, 8, 10 and 11 of the long run, plus velocity vs current from the pulse movies of the same day.

Same sample as Studies 10 (600 mA) and 16 (550 mA), three months earlier. Probably a different region of the sample: the wall here crosses the field of view from top to bottom and moves to −x.

## Data

- `MOBILEDATAMOKE/MOKE 2017/Data2017/2017_02_16/B1743/B1743_Bg1iz500mA2s_1/Pos0`, copied to `~/MOKEdata_local/GdFeCo_B1743_2017_02_16/`.
- 10800 Micro-Manager TIFF frames (`img_000000000_Default_000.tif` ...), 1392×1040, Micro-Manager prefix `TestCronometer_2`.
- **Δt = 2.000 s** from the `ElapsedTime-ms` of the metadata (1.99–2.02 s, no gaps); 6.0 h, starting 16 Feb 2017 21:57.
- **Exposure 20 ms**, against 200 ms in Studies 10 and 16: noisier frames (see `fig_h_mean.png`, right).
- Coil BG1, 500 mA (`Bg1iz500mA`). Field calibration not in the folder.
- Scale 0.1175 µm/px assumed (20x, as Studies 10 and 16; `PixelSizeUm` = 0 in the metadata).
- Camera dark level ~280 counts, domain level ~1050 counts.

## Commands

```
python step2_persistence.py GdFeCo_B1743_500mA_16feb17 --row-mean --eps 0.5 1 1.5 2 3 5 --tau-max 2000
python step7_roughness.py GdFeCo_B1743_500mA_16feb17 --every 36 --rmin 25
python step8_crossover.py GdFeCo_B1743_500mA_16feb17
python step10_universality.py GdFeCo_B1743_500mA_16feb17
python step11_overhangs.py GdFeCo_B1743_500mA_16feb17
```

Step 2 was repeated on the frames where the whole wall moves (`win1500_9500/`, 8000 frames) and on three 360-row pieces of that window (`win1500_9500_rows*/`, `--eps 1.5 3`); each folder holds the cut `h_xt_sub.npy` and its `meta.json`.

## Step 1

`python step1_extract_wall.py ~/MOKEdata_local/GdFeCo_B1743_2017_02_16/B1743_Bg1iz500mA2s_1/Pos0 --outdir GdFeCo_B1743_500mA_16feb17 --um-per-px 0.1175 --rotate 2 --bright-domain --offset auto --flatfield 2 --init-wall --crop-right 1200 --qc-every 200 --no-movie --workers 3`

- **`--rotate 2`:** in the raw frames the bright (growing) domain is on the right and the wall moves to −x. Rotated by 180°, it moves to +x as in Studies 10 and 16. Rows and columns in all outputs refer to the rotated frame.
- **Rows:** 110–1026 (917 rows = 108 µm), the longest GdFeCo wall so far (Study 10: 63 µm, Study 16: 42 µm).
- **Drift:** dy ±0.3 px, dx −1.2…+4.1 px.
- **Motion:** all 10800 frames used. Total advance 57.8 px = 6.8 µm in 6.0 h, **v = 0.0054 px/frame ≈ 0.31 nm/s**. No retreats > 1 px; the warning "mean <h> moves backward in 5080 steps" is sub-pixel noise.
- **Tilt:** 12° → 14°, almost constant (Study 16: 26° → 43°), so step 7 on the base plane should work.
- **Overhangs:** front ≠ area in 4.9 % of the (row, frame) pairs.
- **Kymograph:** avalanche-like steps; the lower rows (y > 800) advance most.
- **Noise:** the frame-to-frame scatter of ⟨h⟩ is larger between frames ~1500 and ~9000 than at the start and end. Not understood yet (illumination?).

**Segmentation check.** The overlays were approved by the user before the analysis.

Figures: `fig_qc_overlay.png`, `fig_kymograph.png`, `fig_h_mean.png`, `fig_tilt_overhang.png`, `fig_drift_arrival.png`; overlays every 1000 frames in `qc/` (the full set, every 200 frames, is local only).

## Velocity vs current: pulse movies of the same day

Eleven short movies in the same folder (`B1743_Bg1iz<I>mA<tp>s_<n>`), taken 16:32–17:47, before the long run, **same field of view and same wall geometry**. Each has 10 frames (5 at 1100 mA), exposure 200 ms, one frame after each field pulse; `tp` in the name is the pulse duration (`Interval_ms` = tp + 8–10 s). That `tp` is the pulse duration is supported by the data: at 800 and 900 mA the displacement per pulse scales with tp.

`pulse_velocity.py` (in this folder) measures the mean displacement u_k after k pulses (fixed threshold at half the domain contrast, swept area / rows crossed) and v = Theil–Sen slope of u_k vs k·tp over pulses 1…n, before the domain reaches the aperture edge. Pulse 1 is left out of the slope (u_1 contains a ~2 µm offset from the threshold even when nothing moves, plus a first-pulse transient). Output `pulse_velocity.txt`, figure `fig_velocity_current.png`.

| I [mA] | tp [s] | frames | u_total [µm] | u_1 [µm] | steps fitted | v [nm/s] | reliable |
|---|---|---|---|---|---|---|---|
| 500 | 90 | 10 | 2.7 | 2.2 | 8 | 0.52 | no |
| 600 | 30 | 10 | 3.9 | 1.9 | 8 | 8.7 | no |
| 700 | 1 | 10 | 2.4 | 2.2 | 8 | 14 | no |
| 700 | 17 | 10 | 13.3 | 3.9 | 8 | 66 | yes |
| 800 | 5 | 10 | 22.1 | 4.8 | 8 | 448 | yes |
| 800 | 10 | 10 | 42.4 | 6.9 | 8 | 439 | yes |
| 900 | 1 | 10 | 23.0 | 4.4 | 8 | 2360 | yes |
| 900 | 2 | 10 | 44.5 | 6.5 | 8 | 2280 | yes |
| 900 | 5 | 10 | 87.3 | 17.8 | 8 | 1730 | yes |
| 1000 | 1 | 10 | 90.0 | 14.4 | 7 | 10900 | yes |
| 1100 | 1 | 5 | 106.1 | 48.0 | 1 | 30900 | no (reaches the edge after 2 pulses) |

- **Creep law over 2.5 decades:** the reliable points (700–1000 mA) give ln v [nm/s] = 62.8 − 301.5 I^(−1/4) [I in mA], a straight line in the creep plot.
- **The steady 500 mA run (this study, 0.31 nm/s) lies on the extrapolation of that line** (predicted 0.38 nm/s), 2.5 decades below the slowest reliable pulse point. Studies 10 (600 mA, 2.4 nm/s) and 16 (550 mA, 0.75 nm/s), measured three months later, lie a factor 2–4 below the line.
- 900 mA with 5 s pulses gives a lower v than 1 and 2 s: the wall gets close to the aperture edge in the last pulses.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted, τ ≤ 2000 frames (67 min), τ in frames of 2 s. σ_Δh = 0.71 px (Study 16: 0.47 px; exposure 20 ms instead of 200 ms).

**Whole run (10800 frames):**

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] |
|---|---|---|---|
| 0.5 | 507 (17 min) | 27 | 0.40 ± 0.02 |
| 1 | 654 (22 min) | 29 | 0.73 ± 0.03 |
| 1.5 | 820 (27 min) | 33 | 1.36 ± 0.05 |
| 2 | 747 (25 min) | 30 | 1.48 ± 0.05 |
| 3 | 808 (27 min) | 32 | 2.10 ± 0.05 |
| 5 | 1347 (45 min) | 42 | 2.10 ± 0.05 |

**The run is not stationary.** At τ = 820, ε = 1.5, Π_t is ≈ 0.55 for start times in the first ~1000 and last ~1800 frames and 0.05–0.1 in between. The white noise of h is the same throughout (0.4 px per point). What changes is how much of the wall moves: in the first ~1500 frames the mean advance is the same as later (≈ 0.35 nm/s), but the median row does not move (the motion is in a few rows, mostly y > 800); after frame ~9500 the wall almost stops (0.15 nm/s). This was the "noisier ⟨h⟩ between frames 1500 and 9000" seen in step 1.

**Moving window, frames 1500–9500** (`win1500_9500/`, v ≈ 0.34 nm/s):

| ε [px] | ε/σ_Δh | τ\* | χ4(τ\*) | ξ(τ\*) fit [µm] | C(n) 1/e length [µm] | ∫C(n)/C(0), n < 6 µm [µm] |
|---|---|---|---|---|---|---|
| 0.5 | 0.7 | 293 (10 min) | 7.4 | 0.21 | | |
| 1 | 1.4 | 319 (11 min) | 12.8 | 0.31 | | |
| 1.5 | 2.1 | 458 (15 min) | 20.8 | 0.98 | 0.59 | 2.6 |
| 2 | 2.7 | 485 (16 min) | 23.7 | 1.38 | 0.82 | 3.0 |
| 3 | 4.1 | 736 (25 min) | 31.0 | 2.03 | 2.11 | 4.1 |
| 5 | 6.9 | 1349 (45 min) | 38.1 | 2.27 | 3.41 | 5.1 |

Same quantities in Studies 16 and 10, for comparison: 1/e length 3.4–3.5 µm (550 mA) and 4.6–4.9 µm (600 mA), integral length 5.5–5.8 and 6.5–6.7 µm, flat in ε from ε/σ_Δh = 3 to 13.

What was checked:
- **χ4 is almost all the uniform term B·L of the fit** (B·L/χ4 ≈ +1.0 at every ε; Study 16: −0.5, Study 10: ≈ 0). C(n) has a positive tail of ≈ 9 % of C(0) out to 35 µm (the largest n computed); the exponential takes only the short part, so ξ(τ\*) from the fit is small.
- **Not image drift.** The drift registration is noisy (σ ≈ 0.35 px per frame in x, exposure 20 ms), but its frame-to-frame steps leak into the median step of h with a slope of 0.02 only.
- **Not a global modulation of the drive** (field or temperature): the advance per 10 min of the four quarters of the wall (27 µm each) is uncorrelated in time (r = −0.24 to 0.37). The wall moves in localized bursts, 10–20 px in 10 min over one quarter while the others rest.
- **Not only the wall length:** on 360-row pieces (the length of Study 16) the B·L term is still positive (B·L/χ4 ≈ +0.5–0.7) and ξ(τ\*) at ε = 1.5 is 1.0, 1.15 and 4.1 µm in the three pieces. χ4 per µm of wall is 0.22–0.31, as in Studies 16 (0.22) and 10 (0.27).
- **Noise matters at small ε:** the correlation lengths grow with ε/σ_Δh here (1/e length 0.6 → 3.4 µm from ε/σ_Δh = 2 to 7), and are flat in Studies 10 and 16. At the largest ε the 1/e length (3.4 µm) is close to the 550 mA value.

**Reading:** ξ(τ\*) from the exponential fit is not comparable with Studies 10 and 16 for this run. Model-free, at ε/σ_Δh ≈ 4–7, the correlation length is 2–3.4 µm (1/e) or 4–5 µm (integral), somewhat below 550 and 600 mA; χ4 per unit length is similar. The positive long-range tail of C(n) (correlated regions of tens of µm) is not explained.

## Roughness and l₀

Step 7 on the base plane (tilt 12–14°, residual tilt after rotation ≤ 1.75°), fit 2.9–27 µm. rms width of the detrended wall 6.2 µm.

| estimator (base plane, fit 2.9–27 µm) | ζ |
|---|---|
| global S(q) | 0.95 ± 0.02 (thirds: 0.99, 0.93, 0.92) |
| local width, rotation per segment | 1.03 ± 0.01 (u_area 1.05; folds excluded 1.06) |
| local S(q, l), l = 12–27 µm | 1.04–1.31 |
| B(r) | 0.81 (biased low) |

- Contour and area give the same ζ to 0.01: few overhangs (front ≠ area in 4.9 % of the pairs, 7 % folded segments at l = 1 µm, none above 12 µm).
- **Crossover l×** (`step8`): **49 µm** from the local width (bootstrap 46–54; variants 0.7–150, 3/30 not bracketed), 43 µm from the global S(q) (41–44), 10 µm from the local slope of w² (8–21). l_tilt = 50 µm, l_fold = 5 µm.
- **Height distribution:** skewness −0.16 ± 0.10 (5 µm windows), −0.19 ± 0.07 (20 µm windows); excess kurtosis +4.1 and +0.6.

## Overhangs (step 11)

P_ov(ℓ) = 0.01 at 5 µm, 0.02 at 11 µm, 0.04 at 23 µm, 0.08 at 41 µm; **P_ov ≈ 0.08–0.10 at l× ≈ 43–49 µm**. P_ov never reaches 1/2 inside the field (l₀_geom > ~100 µm); overhangs > 3 px on the base plane are spaced 446 µm on average. So the 43–49 µm break is **not** l₀ (the overhang scale): it is a candidate for L_av.

## qEW or qKPZ below l×

`step10_universality.py` (`fig_universality.png`), l× = 49 µm.

**A. Roughness:** ζ = 0.97 (0.95–0.98) over [3 µm, l×], 0.94 (0.92–0.95) over [4 µm, l×], 0.91 (0.89–0.92) over [5 µm, l×], folds excluded; contour and u_area agree within 0.02.
- Broken power law over [4 µm, L]: **χ² qEW / qKPZ / equilibrium = 1743 / 6893 / 5775** (unfolded 1774 / 8347 / 6977). Free: ζ₁ = 1.02 (1.01–1.03), with the break pushed to the end of the range. The large χ² of all fixed models (17 points) means none fits well; qEW is the least bad.

**B. KPZ term:**
- Λ from the contour: 0.26–1.66 for L = 2–16 µm (median 0.74), R ≤ 0.18.
- Λ from h(x): 0.39–0.76 (median 0.51), R ≤ 0.09.
- All Λ = O(1), positive; no Λ ≫ 1.

**C. Skewness** −0.1 to −0.2 in 5 and 20 µm windows (Study 10: −0.6).

**Conclusion:** as in Studies 10 and 16, the wall below l× is compatible with qEW (ζ ≈ 0.9–1.05, never near 0.63) with Λ = O(1). The break at 43–49 µm, with P_ov ≈ 0.1 there, is not the overhang scale; it is the first L_av candidate in GdFeCo.

## Caveats

1. **Scale assumed** (0.1175 µm/px) and **field unknown** (only the coil current; H ∝ I assumed in the creep plot).
2. **Exposure 20 ms** in the long run: per-frame noise higher than in Studies 10 and 16.
3. **Pulse velocities:** pulse rise time and the field during the ~8–10 s between pulses are unknown; displacement measured from the swept area with a fixed threshold, not with the step1 pipeline. Points with < 3 µm displacement beyond the first pulse are not reliable.
4. **Different region and day** from Studies 10 and 16: the comparison of v, ξ and ζ between the three runs assumes the same disorder.
5. **Non-stationary motion:** for most of the first ~1500 frames only a few rows move, and the wall almost stops after frame ~9500. Step 2 is reported on both the whole run and the moving window; steps 7–11 use the whole run.
6. **Persistence estimator:** the exponential-plus-constant fit of C(n) does not describe this run (B·L ≈ χ4); ξ(τ\*) is not comparable with Studies 10 and 16 (see above).
7. **l× from one wall** (108 µm, ~2 l×): the break is bracketed by the local width and S(q) estimators, but the fixed-model χ² are all large.
8. **Not run:** steps 3–6 and 9.
