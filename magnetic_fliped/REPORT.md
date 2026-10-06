# Study 1: `magnetic_fliped` (polar MOKE movie, 25 fps)

Data: 78 frames of 768×370 px at 25 fps (Δt = 40 ms), δ ≈ 0.17 µm/px (field of view ≈ 131 × 63 µm).
Bright domain on top, dark domain at the bottom, wall moving upward. The original video
(`magnetic_fliped.mp4`) is not in the repository; `h_xt*.npy` and `meta.json` are, so the analysis
can be rerun from step 2. Pipeline and definitions: see the [README](../README.md).

Commands:

```bash
python step1_extract_wall.py magnetic_fliped.mp4 --um-per-px 0.17
python step2_persistence.py magnetic_fliped --eps 0.25 0.5 0.75 1 1.5 2 3
python step2_persistence.py magnetic_fliped --int
python step3_eps_sweep.py magnetic_fliped
```

## Results

- **Noise:** σ_Δh = 0.25 px, so ε = 3σ_Δh = 0.75 px was chosen (subpixel h). Only 0.02 % of the one-frame displacements are retreats larger than 1 px.
- **τ\*:** χ4(τ) decreases monotonically from τ = 1 for every ε, so **τ\* ≤ 1 frame = 40 ms** and is not resolved at this frame rate. The wall advances ~4 px/frame, so Π(τ) is almost zero by τ ≈ 10.
- **ξ(τ\*)** at ε = 0.75 px: **17.7 ± 0.2 px = 3.01 ± 0.04 µm** (statistical error). The systematic error is ~±3 px (~±0.5 µm), because ξ(1) grows from 14 to 22 px (2.4 to 3.8 µm) as ε goes from 0.25 to 3 px.
- **Mean wall velocity:** 3.95 px/frame = 16.8 µm/s.
- **Shape of ξ(τ):** a plateau from τ = 1 to 2–3 frames, followed by a decay to ~5 px at τ ≈ 7–8.
- **Normalized χ4** (ε = 0.75 px): χ4/[Π(1−Π)] ≈ 45–48 px for τ = 1–3, then decays. It has the same shape as ξ(τ). So the monotonic decay of the raw χ4 comes mostly from Π → 0, not from a loss of cooperativity at short τ. The plateau term B·L accounts for about a third of χ4 at τ = 1.

## Discussion and outlook

### Measuring τ\* needs better time resolution
χ4(τ) is largest at the first measured point, τ = 40 ms, so its peak lies at or below one frame. Seeing χ4 rise, peak and fall requires the first lag to leave most columns persistent, roughly Π(Δt) ≈ 0.8–0.9. Here Π(40 ms) = 0.26, so 74 % of the columns have already moved within a single frame.

As an order-of-magnitude estimate, assume Π(τ) ≈ exp(−τ/τ_p). Then Π(40 ms) = 0.26 gives τ_p ≈ 30 ms, and Π(Δt) ≈ 0.8–0.9 requires Δt ≈ 3–7 ms, that is **~150–300 fps**, 6 to 12 times the current rate. The real decay of Π is probably not a pure exponential, so this is only a guide.

### ξ(τ\*) is already reasonably constrained
ξ(τ) (≈ 17–18 px) and the normalized χ4/[Π(1−Π)] (≈ 45–48 px) both stay on a plateau from τ = 1 to 2–3 frames. If τ\* were far below 40 ms, ξ would already be decaying at τ = 1. The plateau suggests that the maximum of ξ lies near 40–120 ms, so **ξ(τ\*) ≈ 3 µm is a reasonable estimate**. Faster imaging would show whether there is a mild maximum inside the plateau.

The dominant uncertainty in ξ is not temporal but the choice of ε (~±0.5 µm). Reducing it needs a better signal-to-noise ratio in the images, not a higher frame rate.

### Experimental trade-offs
In MOKE imaging, a higher frame rate usually means shorter exposure and noisier images. Noisier images increase σ_Δh, force a larger ε and reduce the spatial resolution of the persistence analysis. Options:

- **Faster imaging with the same signal:** more illumination, a more sensitive camera, or a smaller field of view or binning.
- **Lower driving field:** in the creep regime the wall velocity drops very steeply with field, so avalanches spread over more frames at the same 25 fps. This changes the physical operating point, since ξ and τ\* depend on the field, but measuring ξ(H) and τ\*(H) is itself what connects to creep theory.
- **Pulsed-field protocol:** with short field pulses and one image after each pulse, the effective time resolution is set by the pulse duration rather than by the camera.

### Possible next step: avalanche statistics
Even without resolving τ\*, the current movies allow a direct avalanche analysis. Take Δh between consecutive frames, identify contiguous clusters of columns that moved, and measure the distributions of their sizes (swept area) and lateral extents. With 40 ms between frames some clusters will be several avalanches merged together, but the tails of these distributions are usually still informative. This is not implemented yet; it would be a `step4`.
