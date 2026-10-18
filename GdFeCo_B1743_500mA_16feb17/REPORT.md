# Study 17 (in progress): GdFeCo B1743, 500 mA, 16 Feb 2017

**Status: step 1 done, segmentation awaiting approval; velocity vs current from the pulse movies of the same day.** Steps 2–10 not run yet.

Same sample as Studies 10 (600 mA) and 16 (550 mA), three months earlier. Probably a different region of the sample: the wall here crosses the field of view from top to bottom and moves to −x.

## Data

- `MOBILEDATAMOKE/MOKE 2017/Data2017/2017_02_16/B1743/B1743_Bg1iz500mA2s_1/Pos0`, copied to `~/MOKEdata_local/GdFeCo_B1743_2017_02_16/`.
- 10800 Micro-Manager TIFF frames (`img_000000000_Default_000.tif` ...), 1392×1040, Micro-Manager prefix `TestCronometer_2`.
- **Δt = 2.000 s** from the `ElapsedTime-ms` of the metadata (1.99–2.02 s, no gaps); 6.0 h, starting 16 Feb 2017 21:57.
- **Exposure 20 ms**, against 200 ms in Studies 10 and 16: noisier frames (see `fig_h_mean.png`, right).
- Coil BG1, 500 mA (`Bg1iz500mA`). Field calibration not in the folder.
- Scale 0.1175 µm/px assumed (20x, as Studies 10 and 16; `PixelSizeUm` = 0 in the metadata).
- Camera dark level ~280 counts, domain level ~1050 counts.

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
- **The steady 500 mA run (this study, 0.31 nm/s) lies on the extrapolation of that line** (predicted 0.37 nm/s), 2.5 decades below the slowest reliable pulse point. Studies 10 (600 mA, 2.4 nm/s) and 16 (550 mA, 0.75 nm/s), measured three months later, lie a factor 2–4 below the line.
- 900 mA with 5 s pulses gives a lower v than 1 and 2 s: the wall gets close to the aperture edge in the last pulses.

## Caveats

1. **Scale assumed** (0.1175 µm/px) and **field unknown** (only the coil current; H ∝ I assumed in the creep plot).
2. **Exposure 20 ms** in the long run: per-frame noise higher than in Studies 10 and 16.
3. **Pulse velocities:** pulse rise time and the field during the ~8–10 s between pulses are unknown; displacement measured from the swept area with a fixed threshold, not with the step1 pipeline. Points with < 3 µm displacement beyond the first pulse are not reliable.
4. **Different region and day** from Studies 10 and 16: the comparison of v (and later ξ, ζ) between the three runs assumes the same disorder.
