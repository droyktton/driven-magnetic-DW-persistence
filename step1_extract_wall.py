"""Paso 1: extracción de frames y detección de la pared de dominio h(x,t).

Uso:  python step1_extract_wall.py [video.mp4]
Salidas: frames/, h_xt.npy (entero) y h_xt_sub.npy (subpíxel), n_frames x width en px, fig_qc_overlay.png, fig_h_mean.png
"""
import glob
import os
import subprocess
import sys

import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread
from scipy.ndimage import gaussian_filter1d, median_filter
from skimage.filters import threshold_otsu

VIDEO = sys.argv[1] if len(sys.argv) > 1 else "magnetic_fliped.mp4"
FRAMES_DIR = "frames"
SIGMA_Y = 2.0      # suavizado vertical (px) antes de umbralizar
MEDIAN_X = 5       # kernel de mediana a lo largo de x


def extract_frames():
    os.makedirs(FRAMES_DIR, exist_ok=True)
    for f in glob.glob(os.path.join(FRAMES_DIR, "f_*.png")):
        os.remove(f)
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-v", "error", "-i", VIDEO,
                    "-vsync", "0", os.path.join(FRAMES_DIR, "f_%04d.png")], check=True)


def load_stack():
    files = sorted(glob.glob(os.path.join(FRAMES_DIR, "f_*.png")))
    frames = []
    for f in files:
        im = imread(f).astype(float)
        if im.ndim == 3:
            im = im[..., :3].mean(-1)
        frames.append(im)
    return np.array(frames)          # (T, H, W)


def detect_wall(stack, thr):
    smooth = gaussian_filter1d(stack, SIGMA_Y, axis=1)
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
    h = median_filter(h, size=(1, MEDIAN_X), mode="nearest")
    h_sub = median_filter(h_sub, size=(1, MEDIAN_X), mode="nearest")
    return h, h_sub


def main():
    extract_frames()
    stack = load_stack()
    T, H, W = stack.shape
    thr = threshold_otsu(stack)
    h, h_sub = detect_wall(stack, thr)
    np.save("h_xt.npy", h)
    np.save("h_xt_sub.npy", h_sub)
    print(f"|h_sub - h| medio = {np.abs(h_sub - h).mean():.3f} px")
    print(f"frames={T}  size={W}x{H}  Otsu global={thr:.4f}")
    print(f"columnas sin transición: {(h == H).sum()} de {h.size}")

    # QC: overlay de h(x,t) sobre frames representativos
    idx = np.unique(np.linspace(0, T - 1, 6).round().astype(int))
    fig, axes = plt.subplots(len(idx), 1, figsize=(8, 2.1 * len(idx)))
    for ax, i in zip(axes, idx):
        ax.imshow(stack[i], cmap="gray")
        ax.plot(np.arange(W), h[i], "r-", lw=0.8)
        ax.set_title(f"frame {i}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig("fig_qc_overlay.png", dpi=150)
    plt.close(fig)

    # <h>_x vs frame (posición medida desde abajo para que el avance sea creciente)
    hm = h.mean(axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    axes[0].plot(H - hm, "o-", ms=3)
    axes[0].set_xlabel("frame"); axes[0].set_ylabel(r"$H - \langle h \rangle_x$  [px]")
    axes[0].set_title("Posición media de la pared")
    axes[1].plot(np.arange(1, T), -np.diff(hm), "o-", ms=3)
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_xlabel("frame"); axes[1].set_ylabel(r"$-\Delta\langle h \rangle$ [px/frame]")
    axes[1].set_title("Velocidad media")
    fig.tight_layout()
    fig.savefig("fig_h_mean.png", dpi=150)
    plt.close(fig)

    v = -np.diff(hm)
    print(f"avance total = {hm[0] - hm[-1]:.1f} px; v media = {v.mean():.2f} px/frame; "
          f"pasos negativos (retroceso) = {(v < 0).sum()}")


if __name__ == "__main__":
    main()
