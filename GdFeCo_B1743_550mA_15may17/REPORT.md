# Study 16 (in progress): GdFeCo B1743, 550 mA, 15 May 2017

**Status: step 1 only, segmentation awaiting approval.** Steps 2–10 not run yet.

## Data

- `~/MOKEdata_local/GdFeCo_B1743_2017_05_15/Pos0` (copy of `MOBILEDATAMOKE`), 10000 frames, Interval_ms 2000 (Δt = 2 s, 5.6 h).
- Same sample as Study 10 (600 mA, 16 May 2017), one day earlier.
- Scale 0.1175 µm/px assumed, as in Study 10.

## Step 1

`python step1_extract_wall.py ~/MOKEdata_local/GdFeCo_B1743_2017_05_15/Pos0 --outdir GdFeCo_B1743_550mA_15may17 --um-per-px 0.1175 --bright-domain --offset auto --flatfield 2 --init-wall --crop-right 1200 --qc-every 200 --workers 3`

- Rows with a wall: y = 505–864 (360 rows, 42 µm); above 505 the front runs along the aperture edge.
- Total advance 128.5 px (15 µm) in 5.6 h: **v ≈ 0.75 nm/s** (Study 10, 600 mA: 2.4 nm/s).
- Retreats > 1 px: 0 %. The warning "mean <h> moves backward in 4570 steps" is sub-pixel noise (mean advance 0.013 px/frame).
- Drift < 0.5 px. Front ≠ swept area in 8.6 % of (row, frame) pairs (Study 10: 11 %).
- Tilt 26° → 43°: the lower rows advance much more than the upper ones (Study 10: 27° → 7°). Tilt corrections are not small (1/cos 43° ≈ 1.37).
- Kymograph: avalanche-like steps, tens of rows jumping at once.

Figures: `fig_qc_overlay.png`, `fig_kymograph.png`, `fig_h_mean.png`, `fig_tilt_overhang.png`, `fig_drift_arrival.png`; overlays every 1000 frames in `qc/` (the full set every 200 frames and `overlay.mp4` are local only). `h_xt*.npy` will be committed with the rest of the study.
