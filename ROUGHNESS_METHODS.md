# Local roughness exponent ζ and crossover length l₀: methods, validation and lessons

This is a self-contained description of every method used in this repository to measure the roughness of a domain wall: the local exponent ζ, the crossover length l₀, and the geometric descriptors around them. It covers how each method is defined, what it assumes, what the synthetic tests showed, and where it failed. It is meant as input for testing these estimators on critical configurations of the 2D φ⁴ (Ginzburg–Landau) depinning model of Kolton, Ferrero & Rosso, *Depinning free of the elastic approximation*, [arXiv:2306.13415](https://arxiv.org/abs/2306.13415).

Code: `step7_roughness.py` (global and local estimators, descriptors), `wall_tools.py` (shared helpers), `step8_crossover.py` (l₀), `step9_wavelet.py` (wavelets), `wall_exponents.py` (wavelet and local-frame estimators with their synthetic harness). Results are in `0_8A_15s_20x_RT_1/REPORT.md` (Study 2) and `magnetic_fliped/REPORT.md` (Study 1).

---

## 0. Picture being tested

- **Below l₀:** qEW-like self-affine wall, ζ ≈ 1.25.
- **Above l₀:** ζ_eff ≈ 0.5. Overhangs become important at strong disorder.
- **Scaling of l₀:** l₀ ∼ Δ^−2.2 (Δ: disorder strength).
- **Open question:** at large scales, is the wall self-affine with ζ ≈ 0.5, or a self-similar, rotationally invariant contour? The local tilt spread vs ℓ (§3.1) decides it. In Study 2 it decreases with ℓ (37° → 12° → 3.5°), which favours self-affine.

---

## 1. Preliminaries (common to all methods)

### 1.1 Contour
- **Subpixel contour** of the domain at level 0.5 (`skimage.measure.find_contours`), keeping the longest contour. In the φ⁴ model the natural choice is the φ = 0 level set.
- **Uniform arc length.** The contour is resampled at a uniform arc-length step: 0.5 px for the PCA local width, 1 px for the wavelets.
- **Smooth before thresholding where t must be monotonic.** A contour taken from a binary image is a pixel staircase. Its sideways steps make the coordinate along a window's axis go backwards, so the wavelets discard most windows even for smooth walls. Take the contour of the domain smoothed by σ = 1 px instead (below the optical resolution).

### 1.2 Base plane and heights
- **Base plane.** Rotate the wall by its mean tilt θ₀ so that it is on average parallel to an axis. Rotate the domain as a float image with linear interpolation, threshold at 0.5, then remove a straight line from each frame (the residual tilt). Rotating a binary map with nearest neighbours creates staircase edges and false overhangs.
- **Three single-valued heights per row s** of the rotated frame:
  - **u_area**: the area-conserving column height (left edge + length inside the domain). It is single-valued by construction and equals the true height where there is no overhang. This is the reference.
  - **u_front**: the last crossing.
  - **u_back**: the first crossing.

  Where they differ there is an overhang (gaps ≥ 2 px counted).
- **Defining "the" height of a multivalued wall is ill-posed.** The φ⁴ paper averages the multiple heights, ũ(x); it notes that this creates artificial discontinuities. Methods on segments rotated by their own tilt (§2.2–2.4) avoid choosing a height.

### 1.3 Resolution and imaging
- **The optical PSF smooths the wall and steepens every fluctuation function near the cutoff, which inflates ζ.** Study 2: 0.1175 µm/px, PSF σ ≈ 3–4 px; fits start at rmin = 25 px (2.9 µm). Study 1: 0.17 µm/px, rmin = 16 px.
- **Size of the bias** (synthetic, Gaussian blur σ = 3 px):
  - wavelets: +0.4–0.6 at scales a ≤ 8 px (ℓ ≲ 30 px), still +0.1–0.2 at a = 8–19 px;
  - circle-detrended width: 0.5 → 1.0 for windows of 30–100 points.
- **Recommendation for the model:** run every estimator twice. Once on the ideal φ = 0 contour, and once on a simulated measurement: rasterise at a chosen pixel size, blur with a PSF σ ≈ 3 px, threshold at 0.5, extract the contour as in §1.1.

### 1.4 Errors
- **Time series:** consecutive frames are strongly correlated. Use a block bootstrap over 12 time blocks of frames, refitting ζ for each resample.
- **Along the wall:** cross the time blocks with contiguous segments of the wall (6 per frame → 72 blocks).
- **Systematics dominate:** in every self-test, the spread between variants (fit range, detrending, estimator) was much larger than the bootstrap error. For the model, independent samples replace the bootstrap.

---

## 2. Estimators of ζ

| | estimator | input | output |
|---|---|---|---|
| 2.1 | global S(q), B(r), w(ℓ) | single-valued height on the base plane | one ζ per fit range |
| 2.2 | local width with local rotation (PCA) | contour segments | ⟨w²(ℓ)⟩, ζ_eff(ℓ) |
| 2.3 | local structure factor S(q, ℓ) | contour segments | ζ(ℓ) |
| 2.4 | wavelets DOG-3 / DOG-4 in local frames | contour windows | F(a), ζ_eff(a) |
| 2.5 | local slope of any of the above | — | ζ_eff(ℓ) |

### 2.1 Global estimators on a single-valued height
- **Input:** u(s) (u_area, or h for single-valued walls), with a straight line removed per frame.
- **Structure factor:** S(q) = |FFT(δu · Hann)|² / Σ Hann², averaged over frames, log-binned (8 bins/decade). Fit S ∝ q^−(1+2ζ) over 2π/(L/4) < q < 2π/rmin.
- **Height-difference correlation:** B(r) = ⟨[δu(s+r) − δu(s)]²⟩ ∝ r^{2ζ}, and ζ_eff(r) = ½ d ln B / d ln r. Also computed without pairs near static defects.
- **Line-detrended width:** w²(ℓ), the variance of δu in windows of length ℓ after a linear fit.

**Self-test** (tilted 12°, rotated back):

| true ζ | S(q) | B(r) | w(ℓ) |
|---|---|---|---|
| 0.5 | 0.46 ± 0.08 | 0.44 | 0.52 |
| 0.66 | 0.63 ± 0.16 | 0.56 | 0.64 |
| 1.0 | 0.87 ± 0.14 (0.92 unrotated) | 0.77 | 0.72 |

- **S(q) is the reference among the global estimators.**
- **B(r) and w(ℓ) saturate for ζ ≳ 1** and are useless for the qEW value 1.25.

**Pitfalls**
- **Mixing regimes:** a single fit straddling l₀ gives an exponent of neither regime. Study 2 gave ζ ≈ 0.75–0.82 over 2.9–27 µm while the local estimators show ~1.1 below 10 µm and ~0.5–0.6 above 30 µm.
- **Tilt and overhangs:** without the base plane, the mean tilt dominates the fit (false ζ ≈ 1), and even the sign of the height-distribution skewness flips. Overhangs distort scales up to their size; trust ζ only where u_area, u_front and u_back give the same S(q).

### 2.2 Local width with local rotation (PCA per segment)
- **Segments:** cut each frame's contour into segments of projected length ℓ along the base direction, with windows sliding by ℓ/2. A segment contains every contour point whose projection falls in the window, folds included.
- **Width:** w² is the smallest eigenvalue of the 2×2 covariance of the segment's points. That is the variance normal to the segment's own orthogonal least-squares line, i.e. after rotating the segment by its own tilt.
- **Average:** ⟨w²(ℓ)⟩ over all segments and frames, ∝ ℓ^{2ζ}. Use its local slope (§2.5) rather than a single fit.
- **Variants** (all worth keeping):
  - **tilt filter:** exclude segments whose own tilt differs from the base plane by more than 45°;
  - **fold filter:** exclude segments whose contour length inside the window exceeds 1.5 × the median contour length at the same ℓ. Use the median, not ℓ, as reference, because sub-window wiggles already make contour/ℓ > 1 at every scale (Study 2: 1.17 at 1 µm, 1.52 at 108 µm). A fold-filtered segment is one where the wall turns back on itself, so its best line is meaningless;
  - **controls:** the same PCA on u_area, and the line-detrended width without rotation.
- **Self-test:** 0.53, 0.66, 0.77 for true ζ = 0.5, 0.66, 1.0. Unbiased up to 2/3, **underestimates ζ ≈ 1**.
- **Filters bias by selection:** the filters remove 15–20 % of the segments at small ℓ but none at large ℓ, and the removed ones are wide. Filtering therefore lowers w² only at small ℓ and steepens the curve. Study 2, 3–27 µm: full contour 0.98, folds excluded 1.09, tilt ≤ 45° 1.15.
- **Advantage:** every segment counts and folds need no special treatment.
- **Limitation:** for a strongly bent segment the "rotated width" is a width about one straight line, so bends add to it.

### 2.3 Local structure factor S(q, ℓ)
- **Segments:** the same as in §2.2, each rotated by its own tilt. Project the points onto the segment's own axis and normal, and bin at 1 px along the axis using the mean normal displacement per bin. Folds are averaged; empty bins are interpolated.
- **Non-periodic ends — use end matching, not a window.** A segment does not end at the height it starts, and the FFT sees a jump whose spectrum ∝ q⁻² pulls ζ towards 1/2. Tested on ideal profiles (ζ = 0.5–1.25, ℓ = 100–400 px, fit from the first mode):
  - least-squares line + Hann window: overestimates ζ by 0.05–0.2 in short segments (the window spreads the lowest modes);
  - **subtracting the line through the two end points**, no window: unbiased within ±0.03. The periodic continuation is then continuous, and the remaining slope kink leaks only as q⁻⁴.

  Then S(q, ℓ) = |FFT(profile)|² / N, averaged in log-q bins over segments and frames.
- **Fits:** ζ(ℓ) over 2π/ℓ < q < 2π/rmin (all scales of the segment above the resolution), and over the low band 2π/ℓ < q < 8·2π/ℓ (scales ~ ℓ). Both with all segments and without folded ones. Compare with the global S(q) in the same band.
- **Range:** needs ℓ ≥ 4·rmin. Each fit uses only a few modes (the low band holds modes 1–8), so neighbouring ℓ scatter by ±0.1–0.15, more than the bootstrap error.
- **Full self-test** (tilted synthetic wall → contour → rotated segments): ζ within ~0.1.
- **Artefacts:**
  - a high-q floor that grows with ℓ (probably steps created by the binning where parts of the segment are steep relative to its axis);
  - folded segments make ζ erratic (Study 2: 0.69 to 1.29 at 26–40 µm).

### 2.4 Wavelets DOG-3 / DOG-4 in local frames (`wall_exponents.py`, `step9_wavelet.py`)
- **Curve:** at unit arc-length spacing.
- **Windows:** for each scale a, windows of n = 2H·2.2 ≈ 26a points (H = 6a), each in its own PCA frame (t along the principal axis, u normal to it). Discard windows where t is not strictly monotonic, or that do not cover t_c ± H. Interpolate u(t) on a unit grid.
- **Coefficient:** W_m(a) = Σ ψ_m((t − t_c)/a) u(t) / a, with ψ_m the m-th derivative of a Gaussian: m vanishing moments, L1 normalisation. Then F(a) = √⟨W²⟩ ∼ a^ζ, with ℓ ≈ 4a.
- **Why it is attractive:** the local frame and the vanishing moments make W insensitive to tilt, and to curvature up to order m − 1.
- **Usable scales:**
  - kept fraction of windows ≥ 0.98;
  - window length 26a ≤ 0.6 R_c, where R_c is the 10th percentile of the radius of curvature of the curve smoothed over the window half-length. 0.6 is the largest ratio in the validated "curved" geometry. Where the curve is too short to measure R_c, carry the largest measured R_c forward as a lower bound;
  - 4a ≥ rmin.
- **Fit range:** the longest run of usable scales whose local slopes agree within ±0.15. A second, disjoint plateau is reported if present.
- **Self-test** (12 walls per case, tilted 12°, self-affine below 256 px, slope 6 % at 64 px, contour from the σ = 1 smoothed domain):

| true ζ | no blur, DOG-3 / DOG-4 | blur σ = 3 px, automatic range | blur σ = 3 px, a = 4–8 px |
|---|---|---|---|
| 0.5 | no usable range (too wiggly at the pixel scale) | 0.62 / 0.63 | 1.12 / 1.33 |
| 1.0 | 0.96 ± 0.06 / 0.95 ± 0.07 | 1.12 / 1.12 | 1.53 / 1.72 |
| 1.25 | 1.15 ± 0.11 / 1.14 ± 0.08 | 1.25 / 1.34 | 1.82 / 1.94 |

- **Unbiased without blur:** within ~0.1. The PSF bias is strong.
- **Study 1** (gentle wall, kept ≥ 0.987, window/R_c ≤ 0.54): ζ = 1.09 ± 0.08 (DOG-3) and 1.18 ± 0.07 (DOG-4) over ℓ = 4.6–12.9 µm. This agrees with the local width (1.09) and the global S(q) on the same scales (1.18).
- **Study 2** (tortuous wall): **no usable scale.** The kept fraction is 0.85 at a = 2 px, 0.46 at 8 px and 0.11 at 19 px. The window is longer than 0.6 R_c at every scale. Fitting only the kept windows would select the smoothest parts of the wall, a different fraction at each scale.
- **Negative results (tried, discarded):**
  - **Shorter windows:** kernel support ±4a with the vanishing moments re-imposed by projecting out polynomials of degree < m, windows ~12a points. The kept fraction only rose to 0.94 at a = 2 and 0.42 at 19 px. On synthetic walls tuned to Study 2's kept-fraction curve (tangent angle with 35–45° spread, correlation 15–30 px), ζ came out +0.2 to +0.6 above the same profile on a straight base.
  - **Circular-arc detrending:** Taubin circle fit per window, width from the radial residuals. It needs no monotonic t and keeps all windows, but on the same tortuous walls it is as biased as a straight line (+0.1 to +0.6). The bends sit at scales below the window (15–30 px), so a window holds several bends, not one arc. Arcs only help when the bend radius is much larger than the window, which is the regime where lines already work.
- **Lesson:** on a wall whose bends live at a few µm, those bends are the roughness at those scales. Asking an estimator to "remove" them and return the ζ of an underlying profile is ill-posed. In the φ⁴ model the φ = 0 level set is the object; there is no hidden u.

### 2.5 Local slopes
- **Definition:** ζ_eff(ℓ) = ½ d ln w² / d ln ℓ, or d ln F / d ln a, or −(d ln S / d ln q + 1)/2. Compute it by a linear fit over a window spanning a factor ~3 in scale.
- **Use:** always inspect it before fitting. It shows the plateaus, the crossover and the resolution artefacts (the rise near rmin caused by the PSF).

---

## 3. Geometric descriptors (no exponent, but key for matching)

### 3.1 Local tilt spread vs ℓ
- **Definition:** the standard deviation of the tilt of the segments of §2.2 relative to the base plane, plus the fraction of segments tilted by more than 45°.
- **Self-affine wall with ζ < 1:** the spread decreases with ℓ.
- **Self-similar, rotationally invariant contour:** it stays constant.
- **Study 2:** 37–38° up to ~3 µm, 36° at 3–10 µm, 29° at 10–30 µm, 12° at 30–108 µm, 3.5° at 108 µm. 16–21 % of segments are tilted more than 45° below ~12 µm.
- **Study 1:** 17–19° at a few µm, 7° at 50 µm.

### 3.2 Folded fraction vs ℓ
- **Definition:** the fraction of segments flagged by the fold filter of §2.2.
- **Study 2:** 13–15 % up to ~30 µm, 9 % at 32 µm, 0 above ~50 µm.

### 3.3 Overhang statistics
- **Definition:** on the base plane, the fraction of (row, frame) pairs where u_front ≠ u_back, the extent along the wall and the depth (u_front − u_back) of each overhang.
- **Comparison:** the φ⁴ paper's overhang size u_o grows linearly with L at strong disorder.

### 3.4 Height distribution
- **Definition:** P(δu/σ) of the deviations from the mean position on the base plane of each frame, normalised by the frame's width σ(t). Report skewness and excess kurtosis (block bootstrap), globally and in windows of 5 and 20 µm with a line removed in each.
- **Studies 1 and 2:** the long tail is on the lagging side (Study 2, 20 µm windows: skewness −0.28 ± 0.05).
- **Pitfall:** without the base plane the sign of the skewness can flip.

### 3.5 Wavelet kept fraction and window-to-curvature ratio vs scale
- **Definition:** the two quantities of §2.4, used directly as a measure of tortuosity.
- **Study 2:** kept fraction 0.85 (a = 2 px) → 0.46 (8 px) → 0.11 (19 px), with 26a windows. Window/R_c = 0.85–1.7 at every scale.
- **Study 1:** kept ≥ 0.987, window/R_c ≤ 0.54 up to a = 19 px.

---

## 4. Crossover length l₀ (`step8_crossover.py`)

Three estimates:
- **(a) Local width.** Fit a smooth broken power law to the rotated local width ⟨w²(ℓ)⟩ of §2.2: w² = A ℓ^{2ζ₁}[1 + (ℓ/l₀)^m]^{(2ζ₂−2ζ₁)/m}. The fit is in log space, weighted by the bootstrap error, from ℓ = rmin up.
- **(b) Global S(q).** The same crossover written in q: S = A q^{−(1+2ζ₂)}[1 + (q l₀/2π)^m]^{−(2ζ₁−2ζ₂)/m}, fitted on the log-binned S(q) between 2π/L and 2π/rmin.
- **(c) Model free.** The scale where the local slope ζ_eff(ℓ) of w² first falls below the midpoint between its small-scale value (mean of the first 3 points above rmin) and its large-scale value (mean of the last 3).

**Variants (systematic error).** The main variant is ζ₁ = 1.25, ζ₂ = 0.5 fixed, m = 2, lower cutoff 3 µm. Also: free exponents; m = 1, 2, 4; cutoff 2/3/5 µm; curves full contour / folds excluded / u_area. Statistical error: 200 block-bootstrap refits. If l₀ falls beyond the largest scale, or more than 16 % of the bootstrap fits do, report a lower bound.

**Disorder proxies:**
- l_fold: where the folded fraction falls to half its value at ℓ ≤ 5 µm;
- l_tilt: where the tilt spread falls to half its plateau.

**Self-test.** Walls with S ∝ q^−2 [1 + (q l₀/2π)²]^−0.75 (1.25 → 0.5), l₀ = 100 or 300 px, L = 920 px, 48 walls, with and without blur σ = 3 px. All three estimates recover l₀ within ~±30 %; (c) does worse when l₀ approaches L/3. The bootstrap errors are much smaller than this bias.

**Results.**

| | (a) local width | (b) global S(q) | (c) local slope | proxies |
|---|---|---|---|---|
| Study 2 | 28 µm (bootstrap 25–31; variants 6–47) | 8 µm (6–11; variants 5–44) | 24 µm (8–35) | l_fold 43 µm, l_tilt 39 µm |
| Study 1 | 28 µm (26–29; variants 4–36) | 22 µm (20–25; variants 14–30) | 14 µm (11–22) | l_tilt 31 µm |

- **Both studies:** l₀ ≈ 10–30 µm, with no measurable difference between them.
- **Study 2 estimators disagree by a factor ~3.** The crossover is broad rather than a clean pair of power laws, and the global S(q) mixes the tilts and overhangs that the rotation removes.
- **Free-exponent fits are poorly constrained:** ζ₂ → 0 for S(q).
- **Overhangs and steep segments vanish at ~30–45 µm**, the upper end of l₀.

---

## 5. Summary of ζ by scale (Study 2, rotated estimators)

| ℓ range | local width, full contour | local width, ≤ 45° | local width, folds excluded | local width, u_area | line removed, no rotation | tilt sd |
|---|---|---|---|---|---|---|
| 1–3 µm | 1.12 | 1.40 | 1.33 | 1.30 | 1.09 | 38° |
| 3–10 µm | 1.05 | 1.22 | 1.18 | 1.17 | 0.99 | 36° |
| 10–30 µm | 0.93 | 1.13 | 1.06 | 0.89 | 0.69 | 29° |
| 30–108 µm | 0.57 | 0.60 | 0.69 | 0.38 | 0.57 | 12° |

S(q, ℓ) without folds: 0.9–1.1 for ℓ = 12–40 µm, 0.8–0.9 for 48–72 µm. The global S(q) over the same bands gives 0.77–0.88. Global S(q) over 2.9–27 µm: 0.82 (u_area). Study 1: S(q), w(ℓ), the rotated local width, S(q, ℓ) and the wavelets give ζ ≈ 1.04–1.18 over ~3–33 µm. B(r) gives 0.65, because it saturates. Near the resolution ζ rises to ~1.4–1.6, consistent with PSF smoothing.

---

## 6. Suggested benchmark on φ⁴ critical configurations

1. **Input:**
   - the φ = 0 level set at h_d for several Δ (weak to strong disorder), with many independent samples per Δ. Samples matter more than the number of Δ values;
   - periodic boundary conditions along the wall, if available, so S(q) needs no window.
2. **Two versions of each configuration:**
   - the ideal contour;
   - a simulated measurement: rasterise → PSF σ ≈ 3 px → threshold 0.5 → contour, with the σ = 1 smoothing of §1.1. Choose the pixel size so that l₀ / pixel matches the experiment (l₀ ≈ 85–240 px in Study 2).
3. **For each estimator in §2 and each Δ:**
   - ζ_eff(ℓ) curves;
   - ζ below l₀, to compare with 1.25;
   - l₀ from the three estimates of §4, to compare with l₀ ∼ Δ^−2.2;
   - the effect of the fold and tilt filters;
   - the wavelet kept fraction (at which Δ do the wavelets stop being usable?).
4. **Descriptors of §3 vs Δ:** tilt spread vs ℓ (self-affine or self-similar above l₀?), folded fraction, overhang statistics, height distribution. Compare them with the experimental values listed above to place Studies 1 and 2 on the Δ axis.
5. **Caveat:** the configurations are critical (depinning, no temperature), whereas the experiments are in creep. Creep theory expects depinning-like geometry between the optimal length and the avalanche size, so the comparison is meaningful at those scales, but not identical.
