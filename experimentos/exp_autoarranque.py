import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import linalg as la

# ============================================================================
# AUTO-ARRANQUE: PERTURBACION DE FONDO COMO SEMILLA DEL TIEMPO
#
# La red tiene una carga minima eps0 en cada nodo (energia de vacio).
# gamma ya no es externo: gamma_n = eps0 / (eps0 + E_local)
# Cuando el estado solo tiene energia de fondo, colapsa.
# Cuando tiene mucha energia, evoluciona.
# ============================================================================

print("=" * 70)
print("AUTO-ARRANQUE: PERTURBACION DE FONDO COMO SEMILLA DEL TIEMPO")
print("=" * 70)

N = 32
eps0 = 0.01  # perturbacion de fondo minima

# ------------------------------------------------------------------
# LAPLACIANO
# ------------------------------------------------------------------
L = np.zeros((N, N), dtype=complex)
for i in range(N):
    L[i, i] = 2
    L[i, (i+1)%N] = -1
    L[i, (i-1)%N] = -1

# ------------------------------------------------------------------
# ESTADO INICIAL: PERTURBACION MINIMA
# ------------------------------------------------------------------
# La red no esta perfectamente en reposo.
# Cada nodo vibra con amplitud sqrt(eps0).
psi = np.ones(N, dtype=complex) * np.sqrt(eps0 / N)
psi = psi / np.linalg.norm(psi)

# ------------------------------------------------------------------
# BUCLE PRINCIPAL
# ------------------------------------------------------------------
t_acum = [0.0]
inercia = [0.0]
entropia = [-np.sum(np.abs(psi)**2 * np.log(np.abs(psi)**2 + 1e-15))]
gamma_hist = []
dt_base = 0.01  # paso temporal base minimo

print(f"\nPerturbacion de fondo: eps0 = {eps0}")
print(f"Estado inicial: ruido uniforme de amplitud sqrt(eps0/N) = {np.sqrt(eps0/N):.4f}")
print(f"{'Paso':<6}{'t_acum':<12}{'dt':<12}{'Inercia':<12}{'Entropia':<12}{'gamma_medio':<12}")
print("-" * 66)

for paso in range(1, 200):
    # 1. EVOLUCION UNITARIA
    U = la.expm(-1j * L * dt_base)
    psi = U @ psi
    psi = psi / np.linalg.norm(psi)

    psi_antes = psi.copy()
    P = np.abs(psi) ** 2
    S = -np.sum(P * np.log(P + 1e-15))

    # 2. COLAPSO CON gamma AUTO-AJUSTADO
    # gamma_n depende de la energia local en cada sitio
    E_local = np.abs(psi) ** 2  # densidad de probabilidad ~ energia local
    gamma_n = eps0 / (eps0 + E_local)  # vector de gamma por sitio

    # Sample sitio segun Born
    i_actual = np.random.choice(N, p=P)
    psi_col = np.zeros(N, dtype=complex)
    psi_col[i_actual] = 1.0

    # Colapso: mezcla con peso gamma del sitio actual
    g_efectivo = gamma_n[i_actual]
    psi = (1 - g_efectivo) * psi + g_efectivo * psi_col
    psi = psi / np.linalg.norm(psi)

    # 3. INERCIA
    I = np.sum(np.abs(np.abs(psi_antes)**2 - np.abs(psi)**2))

    # 4. TIEMPO
    t_acum.append(t_acum[-1] + dt_base)
    inercia.append(I)
    entropia.append(S)
    gamma_hist.append(g_efectivo)

    if paso % 30 == 0 or paso == 1:
        print(f"{paso:<6}{t_acum[-1]:<12.4f}{dt_base:<12.6f}{I:<12.6f}{S:<12.4f}{np.mean(gamma_n):<12.4f}")

    if paso == 1:
        print(f"  -> Primer colapso en x={i_actual}, gamma_efectivo = {g_efectivo:.4f}")

    # Criterio de parada
    if I < 1e-10 and paso > 10:
        print(f"\n-> Equilibrio alcanzado en paso {paso}")
        break

# ------------------------------------------------------------------
# ANALISIS
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("ANALISIS: AUTO-ARRANQUE CON DIFERENTES eps0")
print("=" * 70)

def simular_eps0(eps0_val, pasos=200):
    L = np.zeros((N, N), dtype=complex)
    for i in range(N):
        L[i, i] = 2
        L[i, (i+1)%N] = -1
        L[i, (i-1)%N] = -1
    psi = np.ones(N, dtype=complex) * np.sqrt(eps0_val / N)
    psi = psi / np.linalg.norm(psi)
    t = [0.0]; I_h = [0.0]; S_h = [0.0]
    for p in range(1, pasos + 1):
        psi = la.expm(-1j * L * dt_base) @ psi
        psi = psi / np.linalg.norm(psi)
        psi_a = psi.copy()
        P = np.abs(psi) ** 2
        S_h.append(-np.sum(P * np.log(P + 1e-15)))
        E_local = np.abs(psi) ** 2
        gamma_n = eps0_val / (eps0_val + E_local)
        i_act = np.random.choice(N, p=P)
        psi_c = np.zeros(N, dtype=complex); psi_c[i_act] = 1.0
        psi = (1 - gamma_n[i_act]) * psi + gamma_n[i_act] * psi_c
        psi = psi / np.linalg.norm(psi)
        I_h.append(np.sum(np.abs(np.abs(psi_a)**2 - np.abs(psi)**2)))
        t.append(t[-1] + dt_base)
        if I_h[-1] < 1e-10 and p > 10:
            break
    return np.array(t), np.array(I_h), np.array(S_h), p

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

for idx, (ax, eps, color) in enumerate([
    (axes[0,0], 1e-3, '#7c3aff'),
    (axes[0,1], 1e-2, '#00e5a0'),
    (axes[1,0], 1e-1, '#ff4d00'),
    (axes[1,1], 1.0, '#e7298a'),
]):
    t, I, S, p = simular_eps0(eps, 200)
    ax.plot(t, I, color=color, linewidth=2, label=f'eps0={eps}')
    ax.set_xlabel('t', fontsize=10)
    ax.set_ylabel('Inercia', fontsize=10)
    ax.set_title(f'eps0 = {eps}, pasos = {p}', fontsize=11)
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8)
    ax.set_yscale('log')

plt.suptitle('Auto-arranque: gamma emerge de eps0 / (eps0 + E_local)', fontsize=13)
plt.tight_layout()
plt.savefig('C:\\Users\\famil\\Desktop\\serie_armonica\\autoarranque.png', dpi=150)
print("\nGrafico guardado: autoarranque.png")

# ------------------------------------------------------------------
# CONCLUSION
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("CONCLUSION: LA RED SE AUTO-OBSERVA")
print("=" * 70)
print(f"""
Con eps0 = {eps0}:

  - Cada nodo vibra con amplitud minima sqrt(eps0/N) ~ {np.sqrt(eps0/N):.4f}
  - Cuando el estado se extiende uniformemente (solo fondo),
    E_local ~ eps0/N, gamma ~ 1, el sistema colapsa
  - Cuando el estado se concentra en un sitio,
    E_local >> eps0, gamma ~ 0, el sistema evoluciona

  gamma_n = eps0 / (eps0 + |psi(n)|^2)

Esto resuelve el problema de "quien empuja primero":
  La red nunca esta perfectamente en reposo porque eps0 > 0.
  La perturbacion de fondo ES el primer colapso.
  El observador es la red misma comparando su estado local
  con el fondo minimo.

  Ya no hay parametros externos: solo la red, su laplaciano,
  y la carga minima eps0 que la habita.
""")
