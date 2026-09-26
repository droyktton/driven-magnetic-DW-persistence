# Dinámica de una pared de dominio magnético: persistencia y longitud de correlación

Análisis de `magnetic_fliped.mp4`: 78 frames de 768×370 px a 25 fps (Δt = 40 ms). Arriba hay un dominio claro, abajo uno oscuro, y la pared avanza hacia arriba. Escala espacial: δ ≈ 0.17 µm/px, así que el campo de visión mide ≈ 131 × 63 µm.

Requiere numpy, scipy, matplotlib, scikit-image e imageio-ffmpeg. Este último trae el binario de ffmpeg incluido.

## Cómo reproducir

El video original no está en el repositorio. Sin él, se puede arrancar desde el paso 2, porque `h_xt.npy` y `h_xt_sub.npy` sí están incluidos. Para correr el paso 1 hay que copiar `magnetic_fliped.mp4` a esta carpeta.

```bash
python step1_extract_wall.py                       # frames + h(x,t)
python step2_persistence.py --eps 1 2              # h entero
python step2_persistence.py --h h_xt_sub.npy --tag sub_ --eps 0.25 0.5 0.75 1 1.5 2 3
python step3_eps_sweep.py --tag sub_               # sensibilidad a ε
```

La escala es el argumento `--um-per-px` de `step2_persistence.py` y `step3_eps_sweep.py`. Por defecto vale 0.17; con `--um-per-px 0` todo se reporta solo en px.

## Scripts

### `step1_extract_wall.py` — extracción de la pared
1. Extrae todos los frames con `ffmpeg -vsync 0` a `frames/f_XXXX.png`.
2. Carga el stack en escala de grises y calcula un **umbral de Otsu global** sobre todo el stack, para que sea el mismo en todos los frames (resultó 0.459).
3. Suaviza cada columna en vertical (`gaussian_filter1d`, σ = 2 px).
4. Define **h(x,t)** como la primera fila, contando desde arriba, donde la intensidad suavizada cae por debajo del umbral (paso de claro a oscuro).
5. Calcula también una **versión subpíxel**, interpolando linealmente el cruce del umbral entre las filas r−1 y r.
6. Aplica un filtro de mediana a lo largo de x (kernel 5) a ambas versiones para sacar glitches de una sola columna.

Salidas: `h_xt.npy` (entero) y `h_xt_sub.npy` (subpíxel), ambos con forma (n_frames, width) y en px. También genera `fig_qc_overlay.png` y `fig_h_mean.png`.

### `step2_persistence.py` — persistencia, C(n,τ), ξ(τ), χ4(τ)
Toma un archivo h (`--h`) y, para cada ε (`--eps`) y cada τ = 1…n_frames/2, calcula lo siguiente. Todos los promedios son sobre las columnas i y sobre todos los tiempos de inicio t.

- p_i(t,τ) = 1 si |h(x_i,t+τ) − h(x_i,t)| < ε (la columna no se movió), y 0 si se movió.
- Π(τ) = ⟨p_i⟩, la fracción de columnas persistentes.
- C(n,τ) = ⟨p_i p_{i+n}⟩ − Π², calculada con FFT para todo n.
- **ξ(τ)** se obtiene ajustando C(n,τ) = A·e^(−n/ξ) + B. La constante B hace falta porque, al restar el Π² global, la variación de Π entre distintos t deja una meseta en n grande.
- **χ4(τ)** = L·Var_t[Π(t,τ)].
- Como verificación compara tres cosas:
  - χ4, calculado directamente.
  - Σ_n (1−|n|/L)·C(n). Es una identidad exacta, así que coincide siempre.
  - La suma del modelo ajustado, A·coth(1/2ξ) + B·L. Esta es la prueba que no es trivial.
- τ\* es el τ donde χ4 es máximo.

Solo ajusta ξ si hay al menos 200 eventos persistentes (`--min-pers`); si no, deja NaN.

**Estimación de ruido:** el script imprime dos estimadores.
- σ_ruido(x): rugosidad de alta frecuencia a lo largo de x.
- **σ_Δh**: ancho del pico de columnas quietas en Δh(τ = 1). Se mide con los Δh "hacia atrás", que son solo ruido porque la pared nunca retrocede.

Si no se pasa `--eps`, usa ε = 1, 2 y 3σ_Δh.

Con `--um-per-px` (por defecto 0.17), el resumen da ξ y la velocidad media de la pared también en µm y µm/s, y los paneles de ξ tienen un eje derecho en µm.

Salidas por cada ε: `persistence_<tag>eps<ε>.npz`, `fig_Cn_tau_<tag>eps<ε>.png` y `fig_xi_chi4_<tag>eps<ε>.png`. Además imprime una tabla por τ y un resumen.

### `step3_eps_sweep.py` — sensibilidad a ε
Lee todos los `persistence_<tag>eps*.npz` y superpone ξ(τ) y χ4(τ) de cada ε. También grafica ξ(τ=1) y χ4(τ=1) en función de ε. Acepta `--um-per-px` (por defecto 0.17). Salida: `fig_eps_sweep.png`.

## Figuras

| Archivo | Qué muestra |
|---|---|
| `fig_qc_overlay.png` | Control de la segmentación: h(x,t) (entero) en rojo sobre 6 frames repartidos entre el inicio y el final. Sirve para ver que la curva sigue el contorno real. |
| `fig_h_mean.png` | Izquierda: posición media de la pared, H − ⟨h⟩_x, contra el frame, para verificar que el avance es monótono y sin saltos. Derecha: velocidad media −Δ⟨h⟩ en px/frame (~4 px/frame, nunca negativa). |
| `fig_Cn_tau_<tag>eps<ε>.png` | Izquierda: C(n,τ) contra la distancia n, para los primeros τ con ajuste válido. Los puntos son los datos y las líneas el ajuste A·e^(−n/ξ)+B. Derecha: (C−B)/A en escala semilog; si el decaimiento es exponencial se ve una recta de pendiente −1/ξ. |
| `fig_xi_chi4_<tag>eps<ε>.png` | Fila superior en escala lineal, inferior en log-log. Columnas: (1) ξ(τ) con barras de error del ajuste, en px a la izquierda y µm a la derecha; (2) χ4(τ) directo, la suma de C(n) y la suma del modelo, con τ\* marcado en rojo; (3) Π(τ). |
| `fig_eps_sweep.png` | Con h subpíxel. (1) ξ(τ) para cada ε; (2) χ4(τ) en semilog para cada ε; (3) ξ(τ=1) en µm (azul) y χ4(τ=1) (rojo) contra ε. El panel (1) tiene un eje derecho en µm. La franja gris llega hasta ε = 3σ_Δh. |

Nombres de archivo: sin tag, el análisis usa h entero; con `sub_`, usa h subpíxel.

## Resultados principales

- **Ruido:** σ_Δh = 0.25 px, así que se eligió ε = 3σ_Δh = 0.75 px (con h subpíxel).
- **τ\*:** χ4(τ) decrece de forma monótona desde τ = 1 para todos los ε, así que **τ\* ≤ 1 frame = 40 ms** y no se resuelve con este frame rate. La causa es que la pared avanza ~4 px/frame, entonces Π(τ) cae casi a cero en τ ≈ 10.
- **ξ(τ\*)** con ε = 0.75 px: **17.7 ± 0.2 px = 3.01 ± 0.04 µm** de error estadístico, y ~±3 px (~±0.5 µm) de error sistemático porque ξ(1) crece de 14 a 22 px (2.4 a 3.8 µm) cuando ε va de 0.25 a 3 px.
- **Velocidad media de la pared:** 3.95 px/frame = 16.8 µm/s.
- **Forma de ξ(τ):** una meseta entre τ = 1 y 2–3 frames y después una caída hasta ~5 px en τ ≈ 7–8.
