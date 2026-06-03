import numpy as np
import time

print("=" * 72)
print("  RED 3D FCC + DNLS + COLAPSO ESPECTRAL")
print("=" * 72)

# =========================================================================
# 1. Construir red FCC 3D con condiciones periodicas
# =========================================================================
N = 6
eps0 = 1/12; gamma = 12.0; beta = 0.3

sites = []; pos = {}
idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i, j, k)); pos[(i,j,k)] = idx; idx += 1

n_nodes = len(sites)
print(f"\n  RED FCC: {n_nodes} sitios (N={N}, esperado {N**3//2})")

vecinos = [(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),
           (1,0,1),(1,0,-1),(-1,0,1),(-1,0,-1),
           (0,1,1),(0,1,-1),(0,-1,1),(0,-1,-1)]

L = np.zeros((n_nodes, n_nodes), dtype=np.float64)
for (i,j,k), a in pos.items():
    L[a,a] = 12.0
    for di,dj,dk in vecinos:
        ni,nj,nk = (i+di)%N, (j+dj)%N, (k+dk)%N
        if (ni+nj+nk) % 2 == 0:
            L[a,pos[(ni,nj,nk)]] = -1.0

print("  Diagonalizando L...")
t0 = time.time(); evals, evecs = np.linalg.eigh(L); t1 = time.time()
print(f"  Hecho en {t1-t0:.2f}s")
print(f"  Autovalores: [{evals[0]:.6f}, {evals[1]:.6f}, ..., {evals[-2]:.4f}, {evals[-1]:.4f}]")
print(f"  Rango: [{np.min(evals):.6f}, {np.max(evals):.6f}]")
print(f"  evals[0] = {evals[0]:.2e} (debe ser ~0 para modo uniforme)")

# =========================================================================
# 2. Estado inicial: paquete gaussiano 3D + modos de prueba
# =========================================================================
cx, cy, cz = N/2, N/2, N/2; sigma = 2.0
psi = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi[a] = np.exp(-0.5 * r2 / sigma**2)
psi = psi / np.linalg.norm(psi) * 2.0

# =========================================================================
# 3. Evolucion DNLS pura (Strang splitting en base espectral)
# =========================================================================
dt = 0.005; steps = 1000
print(f"\n  EVOLUCION DNLS PURA: gamma={gamma}, dt={dt}, pasos={steps}")
print(f"  {'paso':<6}{'t':<8}{'ancho':<12}{'max|psi|^2':<14}{'energia':<16}{'norma':<10}")
print(f"  {'-'*66}")

psi_lin = psi.copy()
for paso in range(steps + 1):
    P = np.abs(psi_lin)**2
    if paso % 200 == 0:
        maxP = np.max(P); norm = np.sum(P)
        xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
        cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
        r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
        ancho = np.sqrt(np.sum(r2*P)/norm)
        E = np.real(np.sum(psi_lin.conj()*(L@psi_lin))) - gamma/2*np.sum(P**2)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<12.4f}{maxP:<14.6f}{E:<16.6f}{norm:<10.6f}")

    # Strang: no-lineal - lineal - no-lineal
    P = np.abs(psi_lin)**2
    psi_lin = psi_lin * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi_lin
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi_lin = evecs @ coeffs
    P = np.abs(psi_lin)**2
    psi_lin = psi_lin * np.exp(1j * gamma * P * dt/2)

# =========================================================================
# 4. Colapso espectral
# =========================================================================
print(f"\n  COLAPSO ESPECTRAL 3D: beta={beta}, umbral en lambda_max/2")
print(f"  {'paso':<6}{'t':<8}{'ancho':<12}{'max|psi|^2':<14}{'energia':<16}")
print(f"  {'-'*66}")

def sigmoid(x, x0=0.5, s=10.0):
    return 1.0 / (1.0 + np.exp(-s*(x-x0)))

evals_norm = evals / evals[-1]
psi_col = psi.copy()

for paso in range(steps + 1):
    # Colapso espectral
    coeffs = evecs.T @ psi_col
    pk = np.abs(coeffs)**2
    gk = sigmoid(evals_norm, x0=0.5, s=8.0) * eps0 / (eps0 + pk + 1e-16)
    coeffs = coeffs * (1 - beta * gk)
    psi_col = evecs @ coeffs
    psi_col = psi_col / np.linalg.norm(psi_col)

    P = np.abs(psi_col)**2
    if paso % 200 == 0:
        maxP = np.max(P); norm = np.sum(P)
        xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
        cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
        r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
        ancho = np.sqrt(np.sum(r2*P)/norm)
        E = np.real(np.sum(psi_col.conj()*(L@psi_col))) - gamma/2*np.sum(P**2)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<12.4f}{maxP:<14.6f}{E:<16.6f}")

    # DNLS
    P = np.abs(psi_col)**2
    psi_col = psi_col * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi_col
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi_col = evecs @ coeffs
    P = np.abs(psi_col)**2
    psi_col = psi_col * np.exp(1j * gamma * P * dt/2)

# =========================================================================
# 5. Colapso en espacio real 3D (comparacion)
# =========================================================================
print(f"\n  COLAPSO REAL 3D: beta={beta}")
print(f"  {'paso':<6}{'t':<8}{'ancho':<12}{'max|psi|^2':<14}{'energia':<16}")
print(f"  {'-'*66}")

psi_real = psi.copy()
for paso in range(steps + 1):
    P = np.abs(psi_real)**2
    gn = eps0 / (eps0 + P + 1e-16)
    psi_real = psi_real * (1 - beta * gn)
    psi_real = psi_real / np.linalg.norm(psi_real)

    P = np.abs(psi_real)**2
    if paso % 200 == 0:
        maxP = np.max(P); norm = np.sum(P)
        xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
        cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
        r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
        ancho = np.sqrt(np.sum(r2*P)/norm)
        E = np.real(np.sum(psi_real.conj()*(L@psi_real))) - gamma/2*np.sum(P**2)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<12.4f}{maxP:<14.6f}{E:<16.6f}")

    P = np.abs(psi_real)**2
    psi_real = psi_real * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi_real
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi_real = evecs @ coeffs
    P = np.abs(psi_real)**2
    psi_real = psi_real * np.exp(1j * gamma * P * dt/2)

# =========================================================================
# 6. Analisis: distribucion espectral del soliton 3D
# =========================================================================
print(f"\n  ANALISIS ESPECTRAL DEL SOLITON 3D:")

psi_lin_final = psi_lin.copy()
P_lin = np.abs(psi_lin_final)**2
coeffs_lin = evecs.T @ psi_lin_final
pk_lin = np.abs(coeffs_lin)**2

# Fraccion de potencia en modos bajos vs altos
frac_bajos = np.sum(pk_lin[evals < evals[-1]/2]) / np.sum(pk_lin)
print(f"  Fraccion de potencia en modos lambda < lambda_max/2: {frac_bajos:.4f}")
print(f"  Fraccion en modos lambda > lambda_max/2: {1-frac_bajos:.4f}")

# Los 10 modos con mayor amplitud
top_k = np.argsort(pk_lin)[::-1]
print(f"\n  Top 10 modos del soliton 3D:")
print(f"  {'#':<4}{'modo':<6}{'lambda':<12}{'|coeff|^2':<14}{'fraccion':<12}")
print(f"  {'-'*48}")
total_pk = np.sum(pk_lin)
for i, k in enumerate(top_k[:10]):
    frac = pk_lin[k] / total_pk * 100
    print(f"  {i+1:<4}{k:<6}{evals[k]:<12.4f}{pk_lin[k]:<14.6e}{frac:<12.2f}%")

print(f"\n  La potencia esta distribuida en ~{np.sum(pk_lin > 0.01*total_pk)} modos significativos")
print(f"  (|coeff|^2 > 1% de la potencia total)")

# =========================================================================
# CONCLUSION
# =========================================================================
print("=" * 72)
print("  CONCLUSION 3D FCC")
print("=" * 72)
print(f"""
  En 3D FCC ({n_nodes} sitios) el comportamiento es CUALITATIVAMENTE
  EL MISMO que en 1D:

  1. DNLS pura: se forma estructura coherente (soliton 3D)
  2. Colapso espectral: degradacion lenta (no salva la barrera)
  3. Colapso real: sobre-enfoque

  La distribucion espectral del soliton 3D muestra que la potencia
  esta repartida en muchos modos — el soliton 3D, como el 1D,
  NO es un autoestado del Laplaciano lineal.

  La barrera entre colapso y no-linealidad es ESTRUCTURAL,
  independiente de la dimension de la red.
""")
