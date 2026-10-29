# Study 18: GdFeCo B1742 at low temperature (LPS 2019, L. Albornoz)

Data of Lucas Albornoz's PhD, measured at LPS (Orsay) in 2019: GdFeCo film B1742, pMOKE, 20–260 K, pulsed fields
(HV 1–2.5 µs, HP 2 ms, Kepco 2 s). Scale 0.171 µm/px (TEM grid calibration, Nov 2019); H = 0.2115 mT/V for the HV supply.
Local copy: `~/MOKEdata_local/LPS2019_B1742/` (from the LucasHDD disk, `doctorado/resultados experimentales/microscopioPMOKE/2019_atLPS`).
Published analyses of the same film and series:
- Albornoz, Guruciaga, Jeudy, Curiale, Bustingorry, PRB 104, 024203 (2021): roughness;
- Albornoz, Ferrero, Kolton, Jeudy, Bustingorry, Curiale, PRB 104, L060404 (2021): depinning exponents at 10–71 K, which give qEW.

Two kinds of input:
1. **Static walls** (10 series, 20–260 K): the wall pieces u(x) that Lucas extracted for his roughness paper (`*_u.dat`, 17–63 walls per series, 30–260 µm).
2. **Full frames** (6 of the series at 220–260 K): 2048×2048 frames, one per pulse, 36–63 frames per series.

## 1. Roughness of the static walls (`lps_walls.py`, step8, step10 A)

`lps_walls.py` turns the independent walls into step7-style `local_width.npz` and `roughness.npz`, with one bootstrap block per wall. Step10 skips part B for these (`meta.json`: `static_walls`).

| series | walls | l× [µm] (step8 a) | ζ [4 µm, l×] | χ² qEW / qKPZ |
|---|---|---|---|---|
| 20 K, HV 70 V, 1 µs | 51 | 15.2 | 1.09 (1.07–1.11) | 38 / 455 |
| 150 K, HV 245 V, 1 µs | 39 | 16.1 | 1.12 | 15 / 567 |
| 150 K, HV 250 V, 1 µs | 17 | 13.1 | 1.07 | 7 / 312 |
| 220 K, HV 180 V, 2.5 µs | 47 | 22.4 | 1.10 | 15 / 1483 |
| 220 K, HV 200 V, 1 µs | 23 | 21.1 | 1.05 | 9 / 542 |
| 220 K, Kepco 7.06 A, 2 s | 54 | 24.2 | 1.15 | 27 / 1834 |
| 260 K, HP 40 V, 2 ms | 63 | 18.6 | 1.09 | 13 / 2080 |
| 260 K, HV 75 V, 2.5 µs | 59 | 11.3 | 1.02 | 36 / 918 |
| 260 K, HV 85 V, 1 µs | 59 | 10.7 | 0.90 | 171 / 1021 |
| 260 K, Kepco 2.18 A, 2 s | 35 | 27.6 | 1.15 | 7 / 1167 |

- l× is bracketed in all ten series. qEW is favoured over qKPZ and equilibrium in every series, 20 K included. 260 K at 85 V is the weakest: ζ = 0.90, and strict qEW fits badly (χ² 171 against 15 for a free exponent).
- Global S(q) over 3–25 µm: ζ_S = 0.93–1.28, median 1.16.
- Local S(q, ℓ): median 1.10 for ℓ ≤ l×, 0.96 above, 0.91 at ℓ ≥ 3l×.
- The two S(q) results were computed with a one-off script and are not in the repo yet.
- **Published ζ.** The published effective ζ (0.72–0.81) comes from fits of B(r) at small r without local rotation. Our control estimator, without local rotation, is also lower in several series (0.80 at 20 K).
- **Errors.** The walls of one series are probably successive positions of the same wall, so the bootstrap errors are likely underestimated.
- **Hint on l×.** At 260 K the short, strong pulses (HV) give l× ≈ 11 µm, and the long ones at low field (HP 2 ms, Kepco 2 s) give 19–28 µm. At 220 K the Kepco series also gives the largest value. This is the direction expected for L_av, which grows as H drops. However, the HP and Kepco fields are not calibrated here.

## 2. Full frames at 220–260 K (`lps_frames.py`, `lps_frames_post.py`, step11)

**What the frames are.** Each frame is a difference to a reference frame, and shows a bright stripe:
- One edge of the stripe never moves (±0.2–0.4 µm): it is the wall of the reference frame. When the wall crosses it, dark contrast appears on the other side (220 K Kepco, frame 59).
- The other edge is the wall.
- The protocol sweeps the wall back and forth: 4–11 µm per pulse, reversing every 7–17 pulses, so the wall passes 4–11 times through the same mean positions.

**Segmentation (`lps_frames.py`).**
- Gaussian smoothing of 3 px.
- Otsu threshold inside the illuminated disc; the largest bright component is the stripe, with holes filled.
- The frame is rotated so that the stripe runs along x. The two edges are taken per column, refined to subpixel.
- Overlays are in `<series>/frames/qc/`; a few are committed, and the user approved them.
- **220 K HV 200 V is excluded:** the frame is cropped, with no disc and with other domains, and only 3 % of the columns are valid.

**Analysis (`lps_frames_post.py`).**
- The moving edge is the wall.
- **Λ:** step10 B, h(x, t) version, on single pulses (pairs of consecutive frames). The advance is taken with the sign of the pulse. Errors: bootstrap over pulses.
- **Contour for step11:** the stripe contour, minus the reference edge and the field-of-view cuts.
- **Shape memory:** r(ℓ) = ⟨w²_ℓ(h_i − h_j)⟩ / (⟨w²_ℓ(h_i)⟩ + ⟨w²_ℓ(h_j)⟩) for every pair of frames (0 = same shape, 1 = independent shapes).

| series | frames | advance per pulse | Λ, L = 2 / 4 / 8 / 16 µm | max R | P_ov(l×) | l₀^geom | overhang spacing | r(8 / 16 µm), same position, other sweep |
|---|---|---|---|---|---|---|---|---|
| 220 K HV 180 V | 47 | 4.3 µm | 0.52 / 0.52 / 0.73 / 0.81 | 0.11 | 0.04 | > 111 µm | 225 µm | 0.98 / 0.96 |
| 220 K Kepco 2 s | 55 | 5.1 µm | 0.34 / 0.37 / 0.50 / 0.55 | 0.08 | 0.04 | > 96 µm | 424 µm | 0.96 / 0.97 |
| 260 K HP 2 ms | 63 | 5.3 µm | 0.30 / 0.34 / 0.40 / 0.46 | 0.04 | 0.01 | > 111 µm | 5400 µm | 0.86 / 0.81 |
| 260 K Kepco 2 s | 36 | 6.6 µm | 0.20 / 0.17 / 0.13 / 0.23 | 0.03 | 0.01 | > 57 µm | 9000 µm | 0.91 / 0.91 |
| 260 K HV 75 V | 59 | 9.4 µm | 0.16 / 0.15 / 0.16 / 0.36 | 0.01 | 0.00 | > 235 µm | 15500 µm | 0.95 / 0.95 |
| 260 K HV 85 V | 59 | 10.8 µm | 0.27 / 0.28 / 0.33 / 0.42 | 0.02 | 0.00 | > 188 µm | 2100 µm | 0.94 / 0.90 |

Λ bootstrap intervals are typically ±0.05–0.15.

### Persistence and χ₄
Not measurable. In one pulse 97–100 % of the wall moves by more than 0.2 µm: Π(1 pulse) = 0–0.03, and Π = 0 from 2 pulses on. A pulsed field is otherwise equivalent to a continuous one (little happens between pulses), so the limitation is the step per pulse, not the pulses themselves.

### Λ (KPZ term)
- Λ = 0.13–0.81, of order 1 as in Studies 1–17, and R ≤ 0.11.
- Forward and backward pulses give the same values.
- The slowest series (260 K Kepco, 2 s pulses at low field) has the smallest Λ (≈ 0.2), whereas qKPZ would give the largest.
- In 75 and 85 V the step per pulse (9–11 µm) is about l×, so Λ may be biased towards 0 (the wall decorrelates within one pulse). The reliable series are the four with steps of 4–7 µm.
- No qKPZ signature, consistent with the depinning exponents of PRB 104, L060404.

### Overhangs (step11 on the wall)
- P_ov(l×) = 0.00–0.04, and P_ov stays below 1/2 up to the whole wall (57–235 µm).
- The wall is single-valued at the break, so the break at 11–28 µm is not l₀. The candidate is L_av, as in Studies 16 and 17 of GdFeCo.

### Shape memory
- **Consecutive frames of the same sweep** (4–8 µm apart): r drops at large scales, to ≈ 0.5 at 32 µm. The large-scale shape survives one pulse.
- **Same mean position, different sweep:** r ≈ 0.8–1 up to 16 µm. There is almost no memory of the shape: it is not set by the pinning landscape at that position, but by the history of the motion. This is what is expected if the wall visits different metastable states in each sweep.
- At 1–2 µm, r ≈ 1 always, probably because of segmentation noise.
- Pairs are matched by median position only, not column by column.

## Files
- `<series>/`: step8 and step10 A on Lucas's walls.
- `<series>/frames/`:
  - `frames.npz`: h_top and h_bot per frame;
  - `post.npz`, `post.log` and `fig_post.png`: Λ, shape memory, wall sweeps;
  - `overhangs.npz` and `fig_overhangs.png`: step11;
  - `fig_frames.png`;
  - `qc/`: some overlays.
- Not committed (`.gitignore`): `contours.npz` and `contours_stripe.npz`, regenerated by `lps_frames.py` and `lps_frames_post.py`.
