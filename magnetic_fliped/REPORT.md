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

## Roughness

`python step7_roughness.py magnetic_fliped --from-h --every 1 --rmin 16` (`fig_roughness.png`): S(q) of the subpixel h, with a straight line removed from each frame (residual tilt ≤ 4.6°; the wall is close to horizontal, so no rotation is applied). rms width of the detrended wall 2.4 µm. Fitting 2.7–33 µm: **ζ = 1.10 ± 0.04** from S(q) (thirds of the movie: 1.00, 1.21, 1.02; w(ℓ) gives 1.04, B(r) 0.65, which underestimates ζ ≈ 1 according to the self-test of step7). This is much rougher than the slow wall of Study 2 (ζ ≈ 0.75) and closer to the depinning value ζ_dep ≈ 1.25, consistent with a wall driven faster, nearer to depinning; but the two studies differ in sample region, field, frame rate and magnification, so the comparison is only indicative.

**Local width with local rotation** (`fig_local_width.png`): segments of length ℓ rotated by their own tilt give ζ = 1.09 ± 0.02 (line removed without rotation: 1.04). The local tilt is much smaller than in Study 2 (sd 17–19° at a few µm, 7° at 50 µm), so the rotation matters less here and all estimators agree on ζ ≈ 1.05–1.10 (without segments tilted > 45°: 1.10; the fold filter excludes nothing, since h is single-valued in this study).

**Local structure factor S(q, ℓ)** (`fig_local_sq.png`, `local_sq.npz`; same segments, end-point line removed, see the Study 2 report for the method): ζ(ℓ) fitted over all q above the resolution is ≈ 1.0–1.1 for ℓ = 15–88 µm (1.12, 1.07, 0.98, 1.10, 1.08, 1.06, 1.19, 1.02, 1.08, 0.98), consistent with the global ζ = 1.10 ± 0.04; the low band alone (scales ~ ℓ) decreases from ~1.1 at 15–27 µm to 0.76–0.91 at 60–88 µm, as does the global S(q) in the same bands (1.2 → 0.8). The shortest segments (ℓ = 12.6 µm, 2.24) have only four modes above the resolution and are not reliable. No segments are folded here (h is single-valued).

**Height distribution** (`fig_height_distribution.png`; h decreases as the wall advances here, so δh > 0 means lagging behind): on the base plane of each frame, skewness +0.29 ± 0.08 and excess kurtosis −0.42 ± 0.14; in windows of 20 µm, skewness +0.14 ± 0.08, excess kurtosis +1.1 ± 0.4; in windows of 5 µm, excess kurtosis +4.1 ± 1.0. As in Study 2, the long tail is on the lagging side. Without removing the per-frame line (mean only), the skewness changes sign (−0.13): the residual tilt (up to 4.6°) dominates it, which is why the base plane is needed.
**Interpretation.** Following Kolton, Ferrero & Rosso, *Depinning free of the elastic approximation*, [arXiv:2306.13415](https://arxiv.org/abs/2306.13415) (2023), a wall is qEW-like (ζ ≈ 1.2) below a crossover length l₀ that decreases with the disorder strength, with ζ_eff ≈ 0.5 above it. A fit of the crossover (`step8_crossover.py`, `fig_crossover.png`; method in the Study 2 report) gives **l₀ ≈ 28 µm** from the local width (variants 4–36 µm), 22 µm from the global S(q) (14–30 µm) and 14 µm from the local slope; the sd of the local tilt falls to half its plateau at 31 µm. This is the same l₀ as in Study 2 (~10–30 µm) within errors, so the small local tilts of Study 1 do not translate into a measurably larger l₀. h is single-valued by construction in this study, so folds cannot appear. See the Study 2 report for the comparison table.

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
