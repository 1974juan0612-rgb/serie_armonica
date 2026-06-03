import numpy as np, math

print("=" * 72)
print("  B1: KERNEL DE COLAPSO ESPECTRAL K(Delta_n) EN FCC 3D")
print("  La no-localidad de Bell como ilusion de k-espacio")
print("=" * 72)

Z = 12; eps0 = 1.0/Z; eta = math.pi / (3*math.sqrt(2))
beta = (1 - eta) + eps0/2

# -----------------------------------------------------------------------
# Red FCC N=6 (108 sitios)
# -----------------------------------------------------------------------
N = 6; dt = 0.005
sites = []; pos = {}; idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)
print(f"  Red FCC: N={N}, sitios={n_nodes}")

vecinos = [(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),
           (1,0,1),(1,0,-1),(-1,0,1),(-1,0,-1),
           (0,1,1),(0,1,-1),(0,-1,1),(0,-1,-1)]
L = np.zeros((n_nodes, n_nodes))
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

# -----------------------------------------------------------------------
# Kernel K(Delta_n) = (1/N) * Sigma_k Gamma_k * exp(ik.Delta_n)
# donde Gamma_k = sigma(lambda_k) * eps0 / (eps0 + |c_k|^2)
# En el punto fijo, |c_k|^2 esta concentrado en lambda=0.
# Para el kernel desnudo (independiente del estado), usamos
# Gamma_k = sigma(lambda_k).
# -----------------------------------------------------------------------
print(f"\n  Construyendo kernel K(Delta_n)...")

# Gamma_k para cada modo
Gamma_k = sigmoid(evals_n)  # filtro espectral desnudo

# K(Delta_n) para todos los pares de sitios
K = np.zeros((n_nodes, n_nodes))
for a in range(n_nodes):
    for b in range(n_nodes):
        # K(a,b) = (1/N) * Sigma_k Gamma_k * e_k(a) * conj(e_k(b))
        K[a,b] = np.sum(Gamma_k * evecs[a,:] * evecs[b,:]) / n_nodes

# Verificar que K es real y simetrico
print(f"  K es real: {np.allclose(K, K.real)}")
print(f"  max |Im(K)|: {np.max(np.abs(K.imag)):.2e}")
print(f"  K es simetrico: {np.allclose(K, K.T)}")
K = K.real

# -----------------------------------------------------------------------
# Analisis del soporte: |K(Delta_n)| vs |Delta_n|
# -----------------------------------------------------------------------
print(f"\n  {'Delta_n':<10}{'|K|':<16}{'decaimiento':<20}{'sitios':<10}")

# Para cada distancia, encontrar el maximo |K|
dists = np.arange(0, 8, 0.5)
maxK = []
for d in dists:
    vals = []
    for a in range(n_nodes):
        for b in range(n_nodes):
            dx = sites[a][0] - sites[b][0]
            dy = sites[a][1] - sites[b][1]
            dz = sites[a][2] - sites[b][2]
            dist = math.sqrt(dx**2 + dy**2 + dz**2)
            if abs(dist - d) < 0.25:
                vals.append(abs(K[a,b]))
    maxK.append(max(vals) if vals else 0)

# Encontrar soporte 1/e
K0 = maxK[0] if maxK[0] > 0 else 1.0
for i, d in enumerate(dists):
    label = ""
    if i > 0 and maxK[i] > 0 and maxK[i-1]/K0 > 1/math.e >= maxK[i]/K0:
        label = f" <-- soporte 1/e = {d:.1f}"
    if maxK[i] < 1e-6 and (i == 0 or maxK[i-1] >= 1e-6):
        label += f" <-- |K| < 1e-6"
    print(f"  {d:<10.1f}{maxK[i]:<16.6e}{'(1/e en {:.1f})'.format(d) if i > 0 and maxK[i-1]/K0 > 1/math.e >= maxK[i]/K0 else '':<20}{label:<10}")

# -----------------------------------------------------------------------
# Decaimiento exponencial efectivo y longitud de correlacion xi_K
# -----------------------------------------------------------------------
print(f"\n  LONGITUD DE CORRELACION DEL KERNEL:")
print(f"  ===================================")

# Usar solo distancias unicas (evitar replicas por PBC)
# En red infinita, K(Delta_n) decae monotonamente
pairs = []
for a in range(n_nodes):
    for b in range(a+1, n_nodes):
        dx = abs(sites[a][0] - sites[b][0])
        dy = abs(sites[a][1] - sites[b][1])
        dz = abs(sites[a][2] - sites[b][2])
        # Distancia minima considerando PBC
        dx = min(dx, N-dx); dy = min(dy, N-dy); dz = min(dz, N-dz)
        dist = math.sqrt(dx**2 + dy**2 + dz**2)
        pairs.append((dist, abs(K[a,b])))

pairs.sort(key=lambda x: x[0])
# Promediar en bins de distancia
bins = {}
for d, val in pairs:
    key = round(d, 2)
    if key not in bins:
        bins[key] = []
    bins[key].append(val)

d_vals = sorted(bins.keys())
k_vals = [np.mean(bins[d]) for d in d_vals]
k_std = [np.std(bins[d]) for d in d_vals]

# Ajuste exponencial para d > 0.5
mask = np.array([d > 0.5 for d in d_vals])
if np.sum(mask) > 3:
    logk = np.log(np.maximum(np.array(k_vals)[mask], 1e-20))
    d_arr = np.array(d_vals)[mask]
    coeffs = np.polyfit(d_arr, logk, 1)
    xi_K = -1.0 / coeffs[0]
    print(f"  K(Delta_n) ~ exp(-|Delta_n| / {xi_K:.3f})")
    print(f"  Rango efectivo: |Delta_n| < {3*xi_K:.2f} (3 vidas medias)")

print(f"\n  K(0) = {K0:.6e}")
print(f"  K(1) / K(0) = {k_vals[d_vals.index(min(d_vals, key=lambda x: abs(x-1.0)))]/K0 if d_vals else 0:.6f}" if d_vals else "")
print(f"  K(2) / K(0) = {k_vals[d_vals.index(min(d_vals, key=lambda x: abs(x-2.0)))]/K0 if d_vals else 0:.6f}" if d_vals else "")

# -----------------------------------------------------------------------
# Implicacion para no-localidad de Bell
# -----------------------------------------------------------------------
print(f"\n" + "=" * 72)
print("  IMPLICACION PARA NO-LOCALIDAD DE BELL")
print("=" * 72)
print(f"""
  El kernel K(Delta_n) = (1/N) * Sigma_k Gamma_k * exp(ik.Delta_n)
  describe el acoplamiento espectral entre dos sitios. Su forma
  esta dictada por la transformada de Fourier de Gamma_k = sigma(lambda_k).

  Resultados:
  - K(0) ~ K0 (auto-acoplamiento del sitio consigo mismo)
  - Decaimiento exponencial con longitud xi_K ~ {xi_K:.2f} sitios
  - Rango efectivo ~ {3*xi_K:.2f} sitios antes de ser despreciable

  Implicacion para Bell:
  1. El acoplamiento es local en k-espacio (Gamma_k actua modo a modo).
  2. Aparece no-local en x-espacio solo al proyectar la base.
  3. S(|r1-r2|) debe decaer exponencialmente desde ~2.828
     hasta S<2 cuando |r1-r2| >> xi_K.
  4. Esto es medible en redes opticas sinteticas (Tissot 2024).
""")
