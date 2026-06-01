import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import eigsh

# ============================================================================
# EL TIEMPO COMO TASA DE CAMBIO, LA INERCIA COMO OBSERVADOR
# 
# Idea: El tiempo no existe hasta que la incertidumbre se termina.
# - t=0 es el fin de la incertidumbre (colapso a un estado definido)
# - La inercia I = diferencia entre lo probable y lo que ocurre
# - El observador es el que actualiza la probabilidad a realidad
# ============================================================================

print("=" * 70)
print("EL TIEMPO COMO TASA DE CAMBIO CON INERCIA")
print("=" * 70)

# ------------------------------------------------------------------
# 1. RED FCC 1D (LINEA) COMO CASO MAS SIMPLE
# Usamos una cadena lineal con acoplamiento a 1os vecinos
# La "inercia" se define como la discrepancia entre la evolucion
# esperada (por el Laplaciano) y la actual (con colapso)
# ------------------------------------------------------------------
N = 64  # numero de sitios en la red 1D
print(f"\nRed 1D con {N} sitios, acoplamiento a primeros vecinos")

# Laplaciano 1D con condiciones periodicas
diag = 2 * np.ones(N)
off = -1 * np.ones(N - 1)
L = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
L[0, -1] = -1
L[-1, 0] = -1

# ------------------------------------------------------------------
# 2. DEFINICION DEL ESTADO INICIAL
# Antes de t=0: maxima incertidumbre (distribucion uniforme en probabilidad)
# t=0: fin de la incertidumbre -> colapso a un pico en x=0
# ------------------------------------------------------------------
x = np.arange(N)

# Estado inicial ANTES del colapso: gaussiano ancho (incertidumbre maxima)
sigma_0 = N / 4
psi_before = np.exp(-0.5 * ((x - N/2) / sigma_0) ** 2)
psi_before = psi_before / np.sqrt(np.sum(psi_before**2))

# Probabilidad antes del colapso (lo probable)
P_before = psi_before**2

# t=0: colapso -> fin de la incertidumbre
# Despues del colapso: pico en x = N/2
psi_0 = np.zeros(N)
psi_0[N//2] = 1.0  # estado localizado

P_actual = psi_0**2  # lo que ocurre

# La inercia inicial I(0) = diferencia entre lo probable y lo que ocurre
I_0 = np.sum(np.abs(P_before - P_actual))
print(f"Inercia inicial I(0) = {I_0:.6f}")
print(f"(diferencia entre distribucion probable y estado actual tras colapso)")

# ------------------------------------------------------------------
# 3. EVOLUCION TEMPORAL
# Propagamos el estado con la ecuacion de difusion:
#   dpsi/dt = -L psi
# En cada paso, comparamos la evolucion esperada con la actual
# ------------------------------------------------------------------
dt = 0.05
n_steps = 500
t_max = dt * n_steps

psi = psi_0.copy()
P_evol = np.zeros((n_steps + 1, N))
I_hist = np.zeros(n_steps + 1)
entropy_hist = np.zeros(n_steps + 1)
t_hist = np.arange(n_steps + 1) * dt

# t=0
P_evol[0] = psi**2
I_hist[0] = I_0
epsilon = 1e-12
entropy_hist[0] = -np.sum(P_evol[0] * np.log(P_evol[0] + epsilon))

for step in range(1, n_steps + 1):
    # Evolucion con el Laplaciano (difusion en grafo)
    dpsi = -L @ psi
    psi = psi + dt * dpsi
    # Normalizar
    psi = psi / np.sqrt(np.sum(psi**2))

    # Probabilidad actual (lo que ocurre con colapso suave)
    P_actual_step = psi**2

    # Probabilidad esperada (si NO hubiera colapso, solo difusion pura)
    # Es la prediccion del paso anterior: P_evol[step-1]
    P_expected = P_evol[step - 1]

    # Inercia: diferencia entre lo esperado (difusion sin colapso)
    # y lo actual (estado colapsado-normalizado)
    I_step = np.sum(np.abs(P_expected - P_actual_step))

    # Guardar
    P_evol[step] = P_actual_step
    I_hist[step] = I_step
    entropy_hist[step] = -np.sum(P_actual_step * np.log(P_actual_step + epsilon))

print(f"Evolucion completada: {n_steps} pasos, t_max = {t_max:.1f}")

# ------------------------------------------------------------------
# 4. ANALISIS: INERCIA Y EMERGENCIA DEL TIEMPO
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("ANALISIS: INERCIA Y EMERGENCIA DEL TIEMPO")
print("=" * 70)

print(f"""
Definiciones:
  Inercia I(t) = sum_x |P_esperada(x,t) - P_actual(x,t)|
    Donde P_esperada es lo que predice la difusion pura desde t-1
    y P_actual es el estado normalizado tras la evolucion.

  Entropia S(t) = -sum_x P_actual(x,t) log P_actual(x,t)
    Mide la dispersion del estado en la red.

Intuicion:
  - I(0) es maxima: el colapso (fin de incertidumbre) maximiza la inercia
  - I(t) -> 0 cuando el sistema se relaja: lo probable y lo actual convergen
  - El tiempo fluye mientras I(t) > 0. Cuando I(t)=0, el sistema esta en
    equilibrio y el tiempo "se detiene" (estado estacionario).
  - El observador es el proceso que actualiza probabilidad a realidad,
    generando inercia. Sin inercia, no hay direccion temporal.
""")

# ------------------------------------------------------------------
# 5. VISUALIZACION
# ------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# (a) Evolucion de la densidad de probabilidad
ax = axes[0, 0]
im = ax.imshow(P_evol.T, aspect='auto', origin='lower',
               extent=[0, t_max, 0, N], cmap='inferno')
ax.set_xlabel('Tiempo t', fontsize=11)
ax.set_ylabel('Posicion x (sitio de red)', fontsize=11)
ax.set_title('Evolucion de la probabilidad |psi(x,t)|^2', fontsize=12)
plt.colorbar(im, ax=ax)

# (b) Inercia I(t)
ax = axes[0, 1]
ax.plot(t_hist, I_hist, color='#ff4d00', linewidth=2)
ax.axhline(y=0, color='white', linewidth=0.5, alpha=0.3, linestyle='--')
ax.set_xlabel('Tiempo t', fontsize=11)
ax.set_ylabel('Inercia I(t)', fontsize=11)
ax.set_title('Inercia: |P_esperada - P_actual|', fontsize=12)
ax.grid(alpha=0.15)
ax.set_yscale('log')

# (c) Entropia S(t)
ax = axes[1, 0]
ax.plot(t_hist, entropy_hist, color='#00e5a0', linewidth=2)
ax.set_xlabel('Tiempo t', fontsize=11)
ax.set_ylabel('Entropia S(t)', fontsize=11)
ax.set_title('Dispersion del estado en la red', fontsize=12)
ax.grid(alpha=0.15)

# (d) Inercia vs Entropia (diagrama de fase)
ax = axes[1, 1]
ax.plot(entropy_hist, I_hist, color='#7c3aff', linewidth=2)
ax.scatter(entropy_hist[0], I_hist[0], color='#ff4d00', s=80,
           label=f't=0 (colapso)', zorder=5)
ax.scatter(entropy_hist[-1], I_hist[-1], color='#00e5a0', s=80,
           label=f't={t_max:.1f}', zorder=5)
ax.annotate('t=0\nfin incertidumbre',
            (entropy_hist[0], I_hist[0]),
            xytext=(entropy_hist[0]+0.3, I_hist[0]+0.1),
            color='#ff4d00', fontsize=9)
ax.annotate('t -> inf\nequilibrio',
            (entropy_hist[-1], I_hist[-1]),
            xytext=(entropy_hist[-1]+0.3, I_hist[-1]+0.05),
            color='#00e5a0', fontsize=9)
ax.set_xlabel('Entropia S(t)', fontsize=11)
ax.set_ylabel('Inercia I(t)', fontsize=11)
ax.set_title('Diagrama de fase: Inercia vs Entropia', fontsize=12)
ax.legend(fontsize=9)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig('C:\\Users\\famil\\Desktop\\serie_armonica\\tiempo_inercia.png', dpi=150)
print("\nGrafico guardado: tiempo_inercia.png")

# ------------------------------------------------------------------
# 6. CUANDO SE DETIENE EL TIEMPO
# ------------------------------------------------------------------
# El tiempo "se detiene" cuando la inercia es indistinguible de cero.
# En ese punto, la red esta en equilibrio y no hay mas cambio.

umbral = 0.001  # umbral de inercia para considerar "tiempo detenido"
idx_detenido = np.where(I_hist < umbral)[0]

print("\n" + "=" * 70)
print("DETECCION DEL FIN DEL TIEMPO")
print("=" * 70)

if len(idx_detenido) > 0:
    t_detenido = t_hist[idx_detenido[0]]
    print(f"Tiempo hasta equilibrio: t = {t_detenido:.3f}")
    print(f"(Inercia por debajo del umbral {umbral})")
    print(f"Entropia en equilibrio: S = {entropy_hist[idx_detenido[0]]:.4f}")
else:
    print(f"La inercia no baja del umbral {umbral} en t_max = {t_max:.1f}")
    print(f"Inercia final: I({t_max:.1f}) = {I_hist[-1]:.6f}")

# ------------------------------------------------------------------
# 7. RESUMEN CONCEPTUAL
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("RESUMEN: EL TIEMPO Y EL OBSERVADOR")
print("=" * 70)
print(f"""
Sistema: Red lineal de {N} sitios, Laplaciano con difusion
Estado inicial: colapso desde maxima incertidumbre a un pico en x={N//2}

Resultados:
  - Inercia inicial I(0) = {I_0:.4f} (maxima, tras el fin de la incertidumbre)
  - Inercia final I({t_max:.1f}) = {I_hist[-1]:.6f} ({'EN EQUILIBRIO' if I_hist[-1] < umbral else 'AUN EN EVOLUCION'})
  - Entropia inicial S(0) = {entropy_hist[0]:.4f}
  - Entropia final S({t_max:.1f}) = {entropy_hist[-1]:.4f}

Interpretacion conceptual:

  1. t=0 es el fin de la incertidumbre
     -> La probabilidad se actualiza a realidad (colapso)
     -> La inercia I(0) mide el costo de esa actualizacion
  
  2. El tiempo fluye mientras I(t) > 0
     -> Lo probable y lo actual no coinciden
     -> El sistema busca el equilibrio (I -> 0)
  
  3. El observador es la diferencia
     -> No es una entidad externa, es el proceso mismo de
        actualizar probabilidad a realidad en cada paso
     -> Es "posible y necesario" porque sin el no hay
        distincion entre lo que podria ser y lo que es
  
  4. Cuando I(t) = 0, el tiempo se detiene
     -> Estado estacionario de equilibrio
     -> El sistema ha agotado toda posibilidad de cambio
  
  La inercia I(t) = |P_esperada - P_actual| define la
  flecha del tiempo. Es una magnitud medible desde la red.
""")
