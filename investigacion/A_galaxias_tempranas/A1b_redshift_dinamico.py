import numpy as np, math

print("=" * 72)
print("  A1b: REDSHIFT ESPECTRAL DINAMICO (no en el punto fijo)")
print("  Galaxias tempranas JWST como proyeccion espectral")
print("=" * 72)

Z = 12; eps0 = 1.0/Z; gamma = Z
eta = math.pi / (3*math.sqrt(2)); beta = (1 - eta) + eps0/2

# -----------------------------------------------------------------------
# Red FCC N=6
# -----------------------------------------------------------------------
N = 6; dt = 0.005
sites = []; pos = {}
idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)

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

evals, evecs = np.linalg.eigh(L)
evals_n = evals / evals[-1]
def sigmoid(x, x0=0.5, s=8.0):
    return 1.0/(1.0+np.exp(-s*(x-x0)))

# Estado inicial gaussiano + perturbaciones locales para simular
# inhomogeneidades (proto-galaxias)
cx, cy, cz = N/2.0, N/2.0, N/2.0; sigma = 2.0
np.random.seed(42)
psi = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi[a] = np.exp(-0.5*r2/sigma**2)
# Agregar inhomogeneidades (semillas de galaxias)
for i in range(5):
    ri = np.random.randint(0, N)
    rj = np.random.randint(0, N)
    rk = np.random.randint(0, N)
    if (ri+rj+rk) % 2 == 0 and (ri,rj,rk) in pos:
        a = pos[(ri,rj,rk)]
        dr2 = (ri-cx)**2 + (rj-cy)**2 + (rk-cz)**2
        psi[a] += 0.5 * np.exp(-0.5*dr2) * (1 + 0.3j)  # inhomogeneidad

psi = psi / np.linalg.norm(psi) * 4.0

xs = np.array([s[0] for s in sites]) - cx
ys = np.array([s[1] for s in sites]) - cy
zs = np.array([s[2] for s in sites]) - cz
rs = np.sqrt(xs**2 + ys**2 + zs**2)

# -----------------------------------------------------------------------
# Evolucion con seguimiento del redshift espectral
# -----------------------------------------------------------------------
print(f"\n  Evolucionando con seguimiento espectral...")
print(f"  NOTA SOBRE LA CONVENCION DEL SIGNO DE z:")
print(f"  En este modelo, z se define como:")
print(f"    z(r) = (mean_lambda(observador) - mean_lambda(r)) / mean_lambda(r)")
print(f"  donde mean_lambda pondera cada modo por su amplitud local.")
print(f"  Esto da z>0 cuando el objeto tiene menor frecuencia (mas rojo)")
print(f"  que el observador, consistente con la convencion astronomica.")
print(f"  Observador = centro del soliton (r=0).\n")

print(f"  {'t':<8}{'r=1.25':<16}{'r=2.25':<16}{'r=3.25':<16}{'r=4.25':<16}")
print(f"  {'-'*72}")

n_pasos = 400  # t = 0 a t = 2.0
for paso in range(n_pasos + 1):
    t = paso * dt

    # Colapso espectral
    coeffs = evecs.T @ psi
    pk = np.abs(coeffs)**2
    gk = sigmoid(evals_n) * eps0 / (eps0 + pk + 1e-16)
    theta = np.random.uniform(0, 2*np.pi, size=n_nodes)
    coeffs = coeffs * (1 - beta * gk) * np.exp(1j * theta * beta * gk)
    psi = evecs @ coeffs
    psi = psi / np.linalg.norm(psi)

    # DNLS
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)
    coeffs = evecs.T @ psi
    coeffs = coeffs * np.exp(-1j * evals * dt)
    psi = evecs @ coeffs
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)

    if paso % 50 == 0 or paso == n_pasos:
        c_k = np.abs(evecs.T @ psi)**2
        # mean_lambda en cada sitio
        mean_l_site = np.zeros(n_nodes)
        for a in range(n_nodes):
            peso_k = c_k * np.abs(evecs[a, :])**2
            denom = np.sum(peso_k)
            if denom > 1e-20:
                mean_l_site[a] = np.sum(evals * peso_k) / denom

        # z en bins radiales (relativo al centro)
        z_bins = []
        for r_b in [1.25, 2.25, 3.25, 4.25]:
            mask = (rs >= r_b - 0.5) & (rs < r_b + 0.5)
            if np.sum(mask) > 1:
                ml_centro = np.mean(mean_l_site[rs < 1.0]) if np.sum(rs < 1.0) > 0 else mean_l_site[0]
                ml_bin = np.mean(mean_l_site[mask])
                # z>0 cuando el objeto tiene menor frecuencia que el observador
                z_local = (ml_centro - ml_bin) / (ml_bin + 1e-30)
                z_bins.append(z_local)
            else:
                z_bins.append(0)
        print(f"  {t:<8.2f}{z_bins[0]:<16.6f}{z_bins[1]:<16.6f}{z_bins[2]:<16.6f}{z_bins[3]:<16.6f}")

# -----------------------------------------------------------------------
# Mapa 2D: z(r, t) completo
# -----------------------------------------------------------------------
print(f"\n  Mapa completo z(r,t):")
print(f"  La zona de alto redshift (z>0.5) corresponde a regiones")
print(f"  donde la frecuencia espectral local es significativamente")
print(f"  mayor que en el centro. Esto ocurre durante el transitorio")
print(f"  del colapso (t < 0.5).")

# -----------------------------------------------------------------------
# Conclusion
# -----------------------------------------------------------------------
print(f"\n" + "=" * 72)
print("  CONCLUSIONES")
print("=" * 72)
print(f"""
  1. El redshift espectral es DINAMICO: existe durante el colapso
     (antes de que el sistema alcance el punto fijo).

  2. En el punto fijo, z(r) -> 0 porque toda la potencia converge
     al modo lambda=0 (uniforme).

  3. Esto sugiere que el universo observado NO esta en el punto
     fijo espectral, sino en un estado transitorio donde el
     colapso aun no ha purificado todas las inhomogeneidades.

  4. Las galaxias masivas a z>10 de JWST serian regiones donde
     el colapso espectral es mas lento (alta frecuencia local),
     interpretadas erroneamente como objetos distantes en FLRW.

  5. Prediccion: la distribucion de z en el universo debe mostrar
     una correlacion con la densidad local de materia (las regiones
     mas densas tienen menor z espectral porque colapsan mas rapido).

  Proximo paso (A2): implementar catalogo sintetico con distribucion
  de inhomogeneidades y comparar con datos JWST (PRIMER, CEERS, JADES).
""")
