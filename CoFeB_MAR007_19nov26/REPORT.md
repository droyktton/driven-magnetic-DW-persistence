# CoFeB MAR007 (Vincent/Ettore), X20, 19 Nov 2026: slow sequences

CoFeB film, t = 0.7 nm, room temperature. pMOKE creep movies from Ettore, session of 19 Nov (`slow dynamics 11_19`, sequences 37–50). Only the four slow sequences are analysed here. The fast ones (37–45, 21.6–34.3 µT) have not been analysed yet.

## Data

- `~/Downloads/VincentData/transfr-_u_bvtZeq4XL.zip` (35 GB), extracted to `~/Downloads/VincentData/slow dynamics 11_19/`. Logbook: `slow_dynamic_measures_20261119.ods`.
- Camera Hamamatsu ORCA-Flash4.0 (C11440-22CU), 16-bit, 2048×2048. Each saved image is a single **40 ms exposure**; 12–13 Nov used 50 ms. The TIFF tag 270 holds the exposure and the capture time.
- The images are difference images (frame − reference + 32768). Only two reference images exist: `MAR007_ref00037` (used for sequences 37–49) and `MAR007_ref00050`.
- Scale 0.136 µm/px (X20). Field = 0.262 mT/A × coil current.

| seq | current | field | Δt | frames used | rows (wall length) | v |
|---|---|---|---|---|---|---|
| 46 | 0.030 A | 7.9 µT | 2 s | 500 | 937–1783 (115 µm) | 0.023 µm/s |
| 48 | 0.030 A | 7.9 µT | 2 s | 500 | 745–1853 (151 µm) | 0.028 µm/s |
| 49 | 0.030 A | 7.9 µT | 2 s | 499 | 644–1813 (159 µm) | 0.024 µm/s |
| 50 | 0.034 A | 8.9 µT | 2 s | 864 of 1000 | 686–1651 (131 µm) | 0.046 µm/s |

- **seq50:** the logbook says 2000 s at Δt = 4 s, but there are 1000 files 2 s apart (TIFF capture times). It stops at frame 864, where the front reaches the margin.
- **seq47** is a single frame. In the logbook, seq38/39 list 0.02/0.022 V for 0.0905/0.099 A, probably typos for 0.2/0.22 V.

## Segmentation

```
python vincent_walls.py "slow dynamics 11_19" sequenceNN --outdir CoFeB_MAR007_19nov26/seqNN_<I>A \
    --dark-npy DARK --um-per-px 0.136 --fps 0.5 --current <I> --adaptive-thr --stop-frac 0.8 \
    --disk-from "slow dynamics 11_19/MAR007_ref000RR.tif" --core-erode 110
```

- **`--core-erode 110`** (15 µm) keeps the wall away from the border of the field diaphragm (the polygon). With the default 40 px, the line bent near the border.
- **DARK for seq46/48/49:** the mean of the 76 last frames of seq45 (frames 380–455), where the whole disk has switched. seq45 shares the same reference, ref37. Subtracting it removes:
  - the dust and illumination shading that ref37 no longer cancels after the sample drifted;
  - the noise of the single-frame reference.

  Noise per frame drops from 156 to 98 counts. With a constant offset instead (`*_const`, not committed), the segmentation followed the shading.
- **DARK for seq50:** a constant offset of 32768. seq50 has its own reference, ref50.
- The switched domain is on the right and grows to −x, so the default flip applies.

## Analysis

`one.sh` with the same steps as for 12–13 Nov:
- `step2_persistence.py --row-mean --eps 1 2 3 5 8 12`
- `step7_roughness.py --every 5 --rmin 8` (`--every 8` for seq50)
- `step8_crossover.py`
- `step10_universality.py`
- `step11_overhangs.py --tol 12`

The 3 px results of step 11 are kept as `*_tol3`.

| seq | τ* (ε = 5 px) | ξ(τ*) ε = 5 / 8 / 12 px | χ4(τ*) ε = 8 px | ζ local width (PCA) | ζ global S(q) | rms width |
|---|---|---|---|---|---|---|
| 46 | 56 s | 6.0 / 6.4 / 7.1 µm | 15.1 | 0.51 | 0.35 | 1.47 µm |
| 48 | 42 s | 5.5 / 7.0 / 9.2 µm | 13.0 | 0.55 | 0.44 | 1.65 µm |
| 49 | 40 s | 5.6 / 6.5 / 7.7 µm | 15.5 | 0.60 | 0.60 | 1.72 µm |
| 50 | 22 s | 4.3 / 5.2 / 6.3 µm | 12.7 | 0.50 | 0.29 | 1.43 µm |

- **Reproducibility.** The three runs at 7.9 µT agree within about 15 % in v, τ* and ξ. All four sequences lie on the 12 Nov creep line, between seq19 (7.3 µT, 0.017 µm/s) and seq15 (7.9 µT, 0.028 µm/s). So the factor 2–3 offset of 13 Nov and the slow seq23/24 of 12 Nov are well outside the run-to-run scatter.
- **Persistence.** τ* ≈ ε/v to 2ε/v, as on 12–13 Nov. ξ(τ*) is 4–9 µm at ε = 0.7–1.6 µm. At ε = 1–2 px (0.14–0.27 µm), ξ is only 0.2–2.4 µm, which is at the noise level. The ε → 0 limit (L_av) is no better constrained than on 12–13 Nov.
- **Roughness.** ζ is 0.50–0.60 (local width) over 1.1–48 µm, the same as 12–13 Nov. In step 10, qEW has the largest χ² in all four sequences.
- **Overhangs** (tol 1.6 µm): P_ov never reaches 1/2, with one overhang every 1500–1950 µm of wall (none in seq48). The walls are single-valued.
- **Wall-position noise**, from the second differences of h: σ_h = 0.47–0.56 px. On 12–13 Nov it was 0.60–0.63 px (seq16, 19 and 33). The averaged background helps by about 20 % despite 20 % less light. Resolving the scales below L_av needs about 4× less noise. That requires averaging the reference image and accumulating frames over each 2 s interval at acquisition: only 40 ms of every 2 s is used.

## Not committed

- `h_xt.npy`, `h_xt_sub.npy`, `tiff_extra.npz` and `segmentation.gif` (large; regenerate with the command above).
- The superseded segmentations `*_const` (constant background) and `*_erode40` (default margin).
