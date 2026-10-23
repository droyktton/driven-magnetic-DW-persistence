# Study 16: GdFeCo B1743, 550 mA, 15 May 2017

Same sample as Study 10 (600 mA, 16 May 2017), measured one day earlier at a drive 8 % lower. Same pipeline and options as Study 10, except step 7 (see below).

## Data

- `~/MOKEdata_local/GdFeCo_B1743_2017_05_15/Pos0` (copy of `MOBILEDATAMOKE`), 10000 Micro-Manager TIFF frames, Interval_ms 2000 (Δt = 2 s, 5.6 h).
- Scale 0.1175 µm/px assumed, as in Study 10.

## Commands

```
python step1_extract_wall.py ~/MOKEdata_local/GdFeCo_B1743_2017_05_15/Pos0 --outdir GdFeCo_B1743_550mA_15may17 \
    --um-per-px 0.1175 --bright-domain --offset auto --flatfield 2 --init-wall --crop-right 1200 --qc-every 200 --workers 3
python step2_persistence.py GdFeCo_B1743_550mA_15may17 --row-mean --eps 0.5 1 1.5 2 3 5 --tau-max 2000
python step7_roughness.py GdFeCo_B1743_550mA_15may17 --from-h --every 36 --rmin 25
python step8_crossover.py GdFeCo_B1743_550mA_15may17
python step10_universality.py GdFeCo_B1743_550mA_15may17
python step11_overhangs.py GdFeCo_B1743_550mA_15may17   # run 2026-10-23
```

## Wall extraction and checks

**Segmentation check.** `fig_qc_overlay.png` and `qc/overlay_*.png` were approved by the user before the analysis. Overlays are committed every 1000 frames; the full set (every 200 frames) and `overlay.mp4` are local only.

- **Rows:** 505–864 (360 rows = 42 µm). Above row 505 the front runs along the aperture edge. The wall is shorter than in Study 10 (540 rows, 63 µm).
- **Drift:** < 0.5 px.
- **Motion:** all 10000 frames used. Total advance 128.5 px = 15 µm, **mean velocity 0.013 px/frame ≈ 0.75 nm/s**. No retreats > 1 px. Step1 warns that ⟨h⟩ moves backward in 4570 steps; these are sub-pixel (σ_Δh = 0.47 px).
- **Noise:** σ_Δh = 0.47 px, so 3σ_Δh ≈ 1.4 px.
- **Tilt:** θ = 26° at the start, 43° at the end. The lower rows advance much more than the upper ones, so the wall rotates; in Study 10 it straightened (27° → 7°).
- **Overhangs:** the front differs from the area estimate in 8.6 % of the (row, frame) pairs (Study 10: 11 %).
- **Kymograph:** avalanche-like steps, tens of rows jumping at once and then resting.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted, τ ≤ 2000 frames (67 min). τ in frames of 2 s.

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] at τ |
|---|---|---|---|---|
| 0.5 | 431 (14 min) | 6.5 | 3.72 ± 0.11 | 5.01 at 44 min |
| 1 | 454 (15 min) | 7.8 | 4.04 ± 0.09 | 6.08 at 44 min |
| **1.5** | **511 (17 min)** | 9.2 | **4.34 ± 0.12** | 6.88 at 45 min |
| 2 | 499 (17 min) | 9.2 | 4.25 ± 0.11 | 6.80 at 46 min |
| 3 | 524 (17 min) | 9.3 | 4.25 ± 0.11 | 6.59 at 46 min |
| 5 | 679 (23 min) | 9.7 | 4.60 ± 0.14 | 6.13 at 46 min |

Comparison with Study 10 (600 mA):

| | 550 mA (this study) | 600 mA (Study 10) |
|---|---|---|
| v | 0.75 nm/s | 2.4 nm/s |
| τ\* (ε = 1.5) | 17 min | 7.3 min |
| ⟨v⟩τ\* | 6.5 px | 9–13 px |
| ξ(τ\*) (ε ≥ 1) | 4.0–4.6 µm | 4.4–4.5 µm |
| ξ_max | 6.1–6.9 µm | 5.1–5.8 µm |
| χ4(τ\*) | ≈ 9 | ≈ 17 |

- **ξ(τ\*) is the same at both drives (≈ 4.3 µm)**, and almost independent of ε for ε ≥ 1 px.
- **τ\* grows by 2.3 when v drops by 3.2:** the advance per τ\* is somewhat smaller at 550 mA (6.5 px vs 9–13 px).
- **χ4(τ\*) is half that of 600 mA.** Part of this is the shorter wall: χ4 = L Var(Π) counts correlated regions over 360 instead of 540 rows. Per unit length χ4/L is ≈ 0.026 vs 0.032, closer.
- ξ(τ) keeps growing up to τ ≈ 45 min (6–7 µm) and then falls. Beyond τ ≈ 1400 frames, Π < 0.35 and the statistics over the 10000 frames thin out; the fall is not interpreted.

## Roughness and l₀

**Step 7 on the base plane fails here.** With the Study 10 options, step 7 rotates every frame by the mean tilt θ0 = 37°. Because the tilt changes from 26° to 43°, only 82 rows stay inside the rotated image in every frame, with 42 % overhangs on the base plane and no fit range (crash). `step7_rotated_failed.log` keeps that output; its `contours.npz` is the one used by step 10 for Λ.

Step 7 was run with `--from-h` instead: on h(y,t) as measured, without rotation (single-valued, area-based). The local width rotates each segment by its own tilt (PCA), so the global tilt does not enter it; the global S(q) only has a straight line removed per frame (residual tilt up to 47°).

| estimator (`--from-h`, fit 2.9–11 µm) | ζ |
|---|---|
| global S(q) | 0.80 ± 0.12 (thirds: 0.05, 0.70, 1.01) |
| local width, rotation per segment | 0.96 ± 0.02 (segments tilted > 45° excluded: 1.13) |
| local S(q, l), l = 13–30 µm | 0.75–0.85 |
| B(r) | 0.63 (biased low) |

- The first-third S(q) value (0.05) is not physical: early in the run the wall is short in its active part and the straight-line removal over a 26–30° tilted, curved wall dominates. Not investigated further.
- **ζ_eff(ℓ)** of the local width (`fig_crossover.png`, panel 3): ≈ 1.2–1.35 for ℓ = 0.5–3 µm, ≈ 0.9 for 5–15 µm, then falls to 0.1–0.5 above ~20 µm (finite wall, 42 µm). In `fig_universality.png` w²(ℓ) follows the qEW slope (ζ = 1.25) up to ~4 µm.
- **Crossover l₀** (`step8`): 9.2 µm from the local width (bootstrap 7.9–10.5, variants 1.5–150), 7.6 µm from the global S(q) (variants 6.7–60), 26.6 µm from the local slope; l_tilt = 19.1 µm. **l₀ is not determined: ≈ 8–27 µm.**
- **Height distribution** (h, `--from-h`): skewness −0.14 ± 0.12 on the base plane, +0.19 ± 0.10 in 5 µm windows (Study 10: −0.6 ± 0.1); excess kurtosis +3.1 in 5 µm windows.

## Overhangs (step 11, run 2026-10-23)

Run on `contours.npz` of the failed base-plane step 7 (`step7_rotated_failed.log`); step 11 rotates every segment to its own axis, so the global tilt does not enter. P_ov(ℓ) = 0.01 at 5 µm, 0.02 at 11 µm, 0.04 at 22 µm, 0.08 at 39 µm (the largest scale); **P_ov = 0.01–0.02 at the step8 breaks (7.6–9.2 µm)** and never 1/2 (l₀_geom > 39 µm). Overhangs > 3 px on the base plane: one every 600 µm. The wall is single-valued at all measured scales, so whatever bends w²(ℓ) here is not the overhang scale l₀.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A. Roughness below l₀** (l₀ = 9.2 µm): ζ = 0.89 (0.85–0.91) over [3 µm, l₀], 0.83 (0.80–0.86) over [4 µm, l₀], 1.01 (0.95–1.06) over [5 µm, l₀] (3 points). The contour, folds-excluded and u_area columns are identical because step 7 ran on h(y,t); they are not independent checks here.

- Broken power law over [4 µm, L]: **χ² qEW / qKPZ / equilibrium = 44 / 152 / 117**; free: ζ₁ = 0.90 (0.88–0.93) with l₀ = 49 µm, i.e. no break inside the 42 µm wall.

**B. KPZ term:**
- Λ from the contour: 0.27–0.81 for L = 2–8 µm (median 0.46), R ≤ 0.17. The defect mask changes nothing.
- Λ from h(x): 0.4 at 2 µm, 1.2–1.5 at 4 µm, 0.2 at 8 µm (median 0.43).
- All Λ = O(1), within the error bars of isotropic growth; no Λ ≫ 1.

**Conclusion:** compatible with qEW and Λ = O(1), as Study 10. ζ below l₀ is 0.83–0.9 at scales ≥ 3–4 µm, lower than in Study 10 (1.0–1.3); w² follows ζ ≈ 1.25 only below ~4 µm. With the caveats below, this is **not** taken as a drive dependence of ζ.

## Caveats

1. **Scale assumed** (0.1175 µm/px), as in Study 10: µm values, v and l₀ scale with it; ζ, Λ and R do not.
2. **Field unknown:** only the coil current is known; the calibration of coil BG1 is not in the data folder.
3. **Short wall:** 42 µm (360 rows). l₀ cannot be bracketed, and the fit ranges below l₀ hold 3–6 points.
4. **Step 7 without rotation** (`--from-h`): the fixed-angle base plane fails because the tilt grows from 26° to 43°. The global S(q) and B(r) are affected by the tilt and curvature; the local width (rotation per segment) is the reliable estimator. The ζ values are not strictly comparable with Study 10, which used the base plane and the subpixel contour.
5. **Different days and possibly different field of view:** comparing v, τ\* and ζ with Study 10 assumes the same sample region and temperature.
6. **v ratio:** 3.2 for 8 % less current, plausible for creep, but from two runs only.
7. **Not run:** steps 3–6 and 9.

## Related data

- **Velocity vs current of the same sample** (pulse movies of 16 Feb 2017, 500–1100 mA) and the 500 mA steady run: Study 17, [`GdFeCo_B1743_500mA_16feb17/REPORT.md`](../GdFeCo_B1743_500mA_16feb17/REPORT.md) and `fig_velocity_current.png` there. The pulse data follow a creep law, ln v = 62.8 − 301.5 I^(−1/4) (v in nm/s, I in mA); Studies 10 and 16 lie a factor 2–4 below its extrapolation (different day, probably a different region).
