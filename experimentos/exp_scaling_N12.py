import numpy as np
import time

print("=" * 72)
print("  FINITE-SIZE SCALING: N=12 (1728 sitios FCC)")
print("=" * 72)

# =========================================================================
# Construir red FCC N=12
# =========================================================================
N = 12
eps0 = 1/12; gamma = 12.0

sites = []; pos = {}
idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1

n_nodes = len(sites)
print(f"\n  RED FCC: {n_nodes} sitios (esperado {N**3//2})")

vecinos = [(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),
           (1,0,1),(1,0,-1),(-1,0,1),(-1,0,-1),
           (0,1,1),(0,1,-1),(0,-1,1),(0,-1,-1)]

print("  Construyendo L ({0}x{0})...".format(n_nodes))
t0 = time.time()
L = np.zeros((n_nodes, n_nodes), dtype=np.float64)
for (i,j,k), a in pos.items():
    L[a,a] = 12.0
    for di,dj,dk in vecinos:
        ni,nj,nk = (i+di)%N, (j+dj)%N, (k+dk)%N
        if (ni+nj+nk) % 2 == 0:
            L[a,pos[(ni,nj,nk)]] = -1.0
t1 = time.time()
print(f"  Construccion: {t1-t0:.2f}s")

print("  Diagonalizando L...")
t0 = time.time(); evals, evecs = np.linalg.eigh(L); t1 = time.time()
print(f"  Diagonalizacion: {t1-t0:.2f}s")
print(f"  Autovalores: [{evals[0]:.6f}, {evals[1]:.6f}, ..., {evals[-2]:.4f}, {evals[-1]:.4f}]")
print(f"  Rango: [{np.min(evals):.6f}, {np.max(evals):.6f}]")

# =========================================================================
# Espectro: comparar con N=6
# =========================================================================
print(f"\n  ESPECTRO vs N=6:")
print(f"  {'modo':<8}{'N=6 (108)':<14}{'N=12 (864)':<14}{'dif':<10}")
print(f"  {'-'*46}")
# Tomar los primeros 10 y ultimos 5 modos
mostrados = set()
for k in list(range(10)) + list(range(n_nodes-5, n_nodes)):
    if k not in mostrados:
        mostrados.add(k)
        ref = [0, 4, 7, 9, 12, 13, 16]  # autovalores esperados FCC
        if k < len(evals):
            print(f"  {k:<8}{'-':<14}{evals[k]:<14.4f}{'-':<10}")
print(f"\n  Degeneraciones observadas (N=12):")
# Contar degeneraciones
from collections import Counter
evals_round = np.round(evals, 4)
deg = Counter(evals_round)
for val, count in sorted(deg.items()):
    if count >= 4:
        print(f"    lambda = {val:.4f}  x{count}")

# =========================================================================
# Estado inicial y evolucion DNLS pura (sigma=2.0 como N=6)
# =========================================================================
cx, cy, cz = N/2, N/2, N/2; sigma = 2.0
psi = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi[a] = np.exp(-0.5*r2/sigma**2)
psi = psi / np.linalg.norm(psi) * 2.0

dt = 0.005; steps = 500
print(f"\n  DNLS PURA: N=12, gamma={gamma}, dt={dt}, pasos={steps}")
print(f"  {'paso':<6}{'t':<8}{'ancho':<12}{'max|psi|^2':<14}{'energia':<16}")
print(f"  {'-'*66}")

for paso in range(steps + 1):
    if paso % 100 == 0:
        P = np.abs(psi)**2; norm = np.sum(P)
        xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
        cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
        r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
        ancho = np.sqrt(np.sum(r2*P)/norm)
        E = np.real(np.sum(psi.conj()*(L@psi))) - gamma/2*np.sum(P**2)
        maxP = np.max(P)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<12.4f}{maxP:<14.6f}{E:<16.6f}")

    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi = evecs @ coeffs
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)

# =========================================================================
# Colapso espectral N=12
# =========================================================================
beta = 0.3
def sigmoid(x, x0=0.5, s=8.0):
    return 1.0/(1.0+np.exp(-s*(x-x0)))

evals_norm = evals / evals[-1]

# Reiniciar (mismo sigma=2.0)
psi = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi[a] = np.exp(-0.5*r2/sigma**2)
psi = psi / np.linalg.norm(psi) * 2.0

# =========================================================================
# Barrido de amplitud: buscar auto-atrapamiento en N=12
# =========================================================================
for amplitude_norm in [2, 4, 8, 16]:
    psi = np.zeros(n_nodes, dtype=complex)
    for (i,j,k), a in pos.items():
        r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
        psi[a] = np.exp(-0.5*r2/sigma**2)
    psi = psi / np.linalg.norm(psi) * amplitude_norm

    for paso in range(steps + 1):
        if paso % 200 == 0:
            P = np.abs(psi)**2; norm = np.sum(P)
            xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
            cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
            r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
            ancho = np.sqrt(np.sum(r2*P)/norm)
            E = np.real(np.sum(psi.conj()*(L@psi))) - gamma/2*np.sum(P**2)
            maxP = np.max(P)
            if paso == 0:
                print(f"  Norma={amplitude_norm:<4}{'':<10}t=0   ancho={ancho:<8.4f}max={maxP:<10.6f}E={E:<12.4f}")
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        coeffs = evecs.T @ psi
        coeffs = coeffs * np.exp(-1j * evals * dt)
        psi = evecs @ coeffs
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
    P = np.abs(psi)**2; norm = np.sum(P)
    xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
    cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
    r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
    ancho = np.sqrt(np.sum(r2*P)/norm)
    maxP = np.max(P)
    E = np.real(np.sum(psi.conj()*(L@psi))) - gamma/2*np.sum(P**2)
    print(f"  Norma={amplitude_norm:<4}{'':<10}t=2.5 ancho={ancho:<8.4f}max={maxP:<10.6f}E={E:<12.4f}")
    estado = "AUTO-ATRAPADO" if ancho < 4.0 else "EXPANSION"
    print(f"  {'':>16}{'-> ' + estado}")

print(f"\n  COLAPSO ESPECTRAL: N=12, beta={beta}")
print(f"  {'paso':<6}{'t':<8}{'ancho':<12}{'max|psi|^2':<14}{'energia':<16}")
print(f"  {'-'*66}")

psi = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi[a] = np.exp(-0.5*r2/sigma**2)
psi = psi / np.linalg.norm(psi) * 2.0

for paso in range(steps + 1):
    coeffs = evecs.T @ psi
    pk = np.abs(coeffs)**2
    gk = sigmoid(evals_norm, x0=0.5, s=8.0) * eps0 / (eps0 + pk + 1e-16)
    coeffs = coeffs * (1 - beta * gk)
    psi = evecs @ coeffs
    psi = psi / np.linalg.norm(psi)

    if paso % 100 == 0:
        P = np.abs(psi)**2; norm = np.sum(P)
        xs = np.array([s[0] for s in sites]); ys = np.array([s[1] for s in sites]); zs = np.array([s[2] for s in sites])
        cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
        r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
        ancho = np.sqrt(np.sum(r2*P)/norm)
        E = np.real(np.sum(psi.conj()*(L@psi))) - gamma/2*np.sum(P**2)
        maxP = np.max(P)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<12.4f}{maxP:<14.6f}{E:<16.6f}")

    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi = evecs @ coeffs
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)

# =========================================================================
# Analisis espectral final
# =========================================================================
P = np.abs(psi)**2
coeffs = evecs.T @ psi
pk = np.abs(coeffs)**2
total_pk = np.sum(pk)
top_k = np.argsort(pk)[::-1]

frac_bajos = np.sum(pk[evals < evals[-1]/2]) / total_pk
print(f"\n  ANALISIS ESPECTRAL N=12:")
print(f"  Fraccion en modos lambda < lambda_max/2: {frac_bajos:.4f}")

print(f"\n  Top 10 modos del punto fijo:")
print(f"  {'#':<4}{'modo':<6}{'lambda':<12}{'|coeff|^2':<14}{'fraccion':<12}")
print(f"  {'-'*48}")
for i, k in enumerate(top_k[:10]):
    frac = pk[k]/total_pk*100
    print(f"  {i+1:<4}{k:<6}{evals[k]:<12.4f}{pk[k]:<14.6e}{frac:<12.2f}%")

n_signif = np.sum(pk > 0.01*total_pk)
print(f"\n  Modos significativos (>1% potencia): {n_signif}")

# =========================================================================
# Comparacion N=6 vs N=12
# =========================================================================
print("=" * 72)
print("  COMPARACION N=6 vs N=12 (sigma=2.0 en ambos)")
print("=" * 72)
print(f"""
  N=6 (108 sitios):    punto fijo ancho ~2.96, max|psi|^2 ~0.009
  N=12 (864 sitios):   punto fijo ancho ~?, max|psi|^2 ~?

  Si el punto fijo es independiente de N, la barrera es
  ESTRUCTURAL (no un artefacto de tamano finito).
  Si cambia con N, hay que escalar hasta el limite termodinamico.
""")
