"""Paso 1: extracción de frames y detección de la pared de dominio h(x,t).

Uso:  python step1_extract_wall.py VIDEO [--outdir DIR] [--um-per-px 0.17] [--fps F]
                                   [--roi X0 X1 Y0 Y1] [--rotate K] [--flipud] [--fliplr]
                                   [--invert] [--sigma-y 2] [--median-x 5] [--edge-margin 2]

Convención (después de roi/rotate/flip/invert): dominio claro arriba, dominio oscuro
abajo, pared que avanza hacia arriba (h decrece con t). Orden de las transformaciones:
roi (en coordenadas del video original) -> rotate (K x 90° antihorario) -> flipud ->
fliplr -> invert.

Salidas en DIR (por defecto, el nombre del video sin extensión):
  frames/            frames extraídos (ffmpeg -vsync 0)
  h_xt.npy           h(x,t) entero, (n_frames_usados x ancho), en px
  h_xt_sub.npy       h(x,t) subpíxel
  meta.json          parámetros del video y de la detección (los leen step2 y step3)
  fig_qc_overlay.png, fig_h_mean.png
"""
import argparse
import glob
import json
import os
import subprocess

import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread
from scipy.ndimage import gaussian_filter1d, median_filter
from skimage.filters import threshold_otsu

ap = argparse.ArgumentParser()
ap.add_argument("video")
ap.add_argument("--outdir", default=None, help="carpeta de salida (default: nombre del video)")
ap.add_argument("--um-per-px", type=float, default=0.17,
                help="escala espacial en µm/px, se guarda en meta.json (0 = desconocida)")
ap.add_argument("--fps", type=float, default=None, help="frame rate real (default: el del video)")
ap.add_argument("--roi", type=int, nargs=4, metavar=("X0", "X1", "Y0", "Y1"), default=None,
                help="recorte en px del video original: columnas X0:X1, filas Y0:Y1")
ap.add_argument("--rotate", type=int, default=0, help="rotar K x 90° antihorario")
ap.add_argument("--flipud", action="store_true", help="invertir arriba/abajo")
ap.add_argument("--fliplr", action="store_true", help="invertir izquierda/derecha")
ap.add_argument("--invert", action="store_true", help="invertir el contraste (claro <-> oscuro)")
ap.add_argument("--sigma-y", type=float, default=2.0, help="suavizado vertical (px) antes de umbralizar")
ap.add_argument("--median-x", type=int, default=5, help="kernel de mediana a lo largo de x")
ap.add_argument("--edge-margin", type=int, default=2,
                help="un frame es válido solo si en todas las columnas la pared está a más "
                     "de este número de px de los bordes superior e inferior")
ap.add_argument("--retreat-px", type=float, default=1.0,
                help="un Delta h(tau=1) mayor que esto cuenta como retroceso real")
args = ap.parse_args()

outdir = args.outdir or os.path.splitext(os.path.basename(args.video))[0]
frames_dir = os.path.join(outdir, "frames")


def out(name):
    return os.path.join(outdir, name)


def video_fps(path):
    gen = imageio_ffmpeg.read_frames(path)
    meta = next(gen)
    gen.close()
    return float(meta["fps"])


def extract_frames():
    os.makedirs(frames_dir, exist_ok=True)
    for f in glob.glob(os.path.join(frames_dir, "f_*.png")):
        os.remove(f)
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-v", "error", "-i", args.video,
                    "-vsync", "0", os.path.join(frames_dir, "f_%04d.png")], check=True)


def load_stack():
    files = sorted(glob.glob(os.path.join(frames_dir, "f_*.png")))
    frames = []
    for f in files:
        im = imread(f).astype(float)
        if im.ndim == 3:
            im = im[..., :3].mean(-1)
        frames.append(im)
    stack = np.array(frames)          # (T, H, W)
    if args.roi:
        x0, x1, y0, y1 = args.roi
        stack = stack[:, y0:y1, x0:x1]
    if args.rotate % 4:
        stack = np.rot90(stack, args.rotate, axes=(1, 2))
    if args.flipud:
        stack = stack[:, ::-1, :]
    if args.fliplr:
        stack = stack[:, :, ::-1]
    if args.invert:
        stack = stack.max() + stack.min() - stack
    return np.ascontiguousarray(stack)


def detect_wall(stack, thr):
    smooth = gaussian_filter1d(stack, args.sigma_y, axis=1)
    below = smooth < thr
    H = stack.shape[1]
    # primera fila (desde arriba) por debajo del umbral; H si la columna es toda clara
    r = np.where(below.any(axis=1), below.argmax(axis=1), H)
    h = r.astype(float)
    # subpíxel: interpolación lineal del cruce del umbral entre las filas r-1 y r
    ok = (r > 0) & (r < H)
    t_idx, x_idx = np.nonzero(ok)
    rr = r[ok]
    s0 = smooth[t_idx, rr - 1, x_idx]
    s1 = smooth[t_idx, rr, x_idx]
    h_sub = h.copy()
    h_sub[ok] = rr - 1 + (s0 - thr) / (s0 - s1)
    h = median_filter(h, size=(1, args.median_x), mode="nearest")
    h_sub = median_filter(h_sub, size=(1, args.median_x), mode="nearest")
    return r, h, h_sub


def ranges(idx):
    """[0,1,2,5,7,8] -> '0-2, 5, 7-8'"""
    out, start = [], None
    for i, v in enumerate(idx):
        if start is None:
            start = v
        if i + 1 == len(idx) or idx[i + 1] != v + 1:
            out.append(f"{start}" if start == v else f"{start}-{v}")
            start = None
    return ", ".join(out)


def longest_run(mask):
    """(inicio, fin) del tramo contiguo más largo con mask == True (fin exclusivo)."""
    best, start = (0, 0), None
    for i, m in enumerate(list(mask) + [False]):
        if m and start is None:
            start = i
        elif not m and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return best


def main():
    os.makedirs(outdir, exist_ok=True)
    fps = args.fps or video_fps(args.video)
    extract_frames()
    stack = load_stack()
    T, H, W = stack.shape
    thr = threshold_otsu(stack)
    r, h, h_sub = detect_wall(stack, thr)
    print(f"video={args.video}  frames={T}  tamaño (tras roi/rotación)={W}x{H}  "
          f"fps={fps:g}  Otsu global={thr:.4f}")

    # --- chequeo de bordes: la pared tiene que estar dentro del cuadro en todas las columnas.
    # Si toca el borde superior (o no aparece), h queda clavada y esas columnas parecerían
    # "quietas", inflando Pi y chi4. Se usa el tramo contiguo más largo de frames válidos.
    m = args.edge_margin
    frame_ok = ((r >= m) & (r < H - m)).all(axis=1)
    t0, t1 = longest_run(frame_ok)
    warnings = []
    if t1 - t0 < T:
        bad = np.where(~frame_ok)[0]
        warnings.append(f"{T - (t1 - t0)} frames descartados: la pared toca el borde o no "
                        f"aparece en frames {ranges(bad.tolist())}; se usan los frames {t0}..{t1 - 1}")
    if t1 - t0 < 10:
        raise SystemExit(f"ERROR: solo {t1 - t0} frames válidos consecutivos; revisá la "
                         f"orientación (--rotate/--flipud/--invert) o --roi.")
    h, h_sub = h[t0:t1], h_sub[t0:t1]

    # --- chequeo de avance monótono (el estimador de ruido de step2 lo supone)
    hm = h_sub.mean(axis=1)
    v = -np.diff(hm)
    dh1 = h_sub[1:] - h_sub[:-1]
    frac_retreat = float((dh1 > args.retreat_px).mean())
    if hm[0] - hm[-1] <= 0:
        warnings.append("la pared NO avanza hacia arriba en promedio: revisá la orientación "
                        "(--flipud/--rotate) o el contraste (--invert)")
    if (v < 0).any():
        warnings.append(f"el promedio <h> retrocede en {(v < 0).sum()} pasos")
    if frac_retreat > 0.01:
        warnings.append(f"{100 * frac_retreat:.1f}% de los Delta h(tau=1) son retrocesos "
                        f"> {args.retreat_px:g} px: el eps automático de step2 (3 sigma_Dh) "
                        f"no es confiable; pasá --eps a mano")

    np.save(out("h_xt.npy"), h)
    np.save(out("h_xt_sub.npy"), h_sub)
    meta = dict(video=os.path.abspath(args.video), fps=fps,
                um_per_px=args.um_per_px or None, n_frames_video=int(T),
                frame_first=int(t0), frame_last=int(t1 - 1), height=int(H), width=int(W),
                otsu_threshold=float(thr), roi=args.roi, rotate=args.rotate,
                flipud=args.flipud, fliplr=args.fliplr, invert=args.invert,
                sigma_y=args.sigma_y, median_x=args.median_x, edge_margin=m,
                frac_retreat=frac_retreat, warnings=warnings)
    with open(out("meta.json"), "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    # QC: overlay de h(x,t) sobre frames representativos del tramo usado
    idx = np.unique(np.linspace(t0, t1 - 1, 6).round().astype(int))
    fig, axes = plt.subplots(len(idx), 1, figsize=(8, 2.1 * len(idx) * max(H / W, 0.3) / 0.48))
    for ax, i in zip(np.atleast_1d(axes), idx):
        ax.imshow(stack[i], cmap="gray")
        ax.plot(np.arange(W), h_sub[i - t0], "r-", lw=0.8)
        ax.set_title(f"frame {i}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out("fig_qc_overlay.png"), dpi=150)
    plt.close(fig)

    # <h>_x vs frame (posición medida desde abajo para que el avance sea creciente)
    frames = np.arange(t0, t1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    axes[0].plot(frames, H - hm, "o-", ms=3)
    axes[0].set_xlabel("frame"); axes[0].set_ylabel(r"$H - \langle h \rangle_x$  [px]")
    axes[0].set_title("Posición media de la pared")
    axes[1].plot(frames[1:], v, "o-", ms=3)
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_xlabel("frame"); axes[1].set_ylabel(r"$-\Delta\langle h \rangle$ [px/frame]")
    axes[1].set_title("Velocidad media")
    for ax in axes:
        ax.set_xlim(-0.5, T - 0.5)
        for a, b in ((0, t0), (t1, T)):   # frames descartados
            if b > a:
                ax.axvspan(a - 0.5, b - 0.5, color="0.85", zorder=0)
    fig.tight_layout()
    fig.savefig(out("fig_h_mean.png"), dpi=150)
    plt.close(fig)

    print(f"frames usados: {t0}..{t1 - 1} ({t1 - t0} de {T})")
    print(f"avance total = {hm[0] - hm[-1]:.1f} px; v media = {v.mean():.2f} px/frame; "
          f"retrocesos > {args.retreat_px:g} px: {100 * frac_retreat:.2f}%")
    for w in warnings:
        print("ADVERTENCIA:", w)
    print(f"salidas en {outdir}/")


if __name__ == "__main__":
    main()
