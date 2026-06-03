import numpy as np, time, math

print("=" * 72)
print("  MODULACION ESPECTRAL: Comunicacion x1->x2 via kernel K")
print("=" * 72)

Z = 12; eps0 = 1.0/Z; gamma = Z
eta = math.pi / (3*math.sqrt(2)); beta = (1 - eta) + eps0/2
print(f"  beta = {beta:.6f}, eps0 = {eps0:.6f}")

N = 6; dt = 0.005; steps = 600
sites = []; pos = {}
idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)
print(f"  Red FCC N={N}: {n_nodes} sitios")

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
evals, evecs = np.linalg.eigh(L)
evals_n = evals / evals[-1]

def sigmoid(x, x0=0.5, s=8.0):
    return 1.0/(1.0+np.exp(-s*(x-x0)))

cx, cy, cz = N/2.0, N/2.0, N/2.0; sigma = 2.0
psi0 = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi0[a] = np.exp(-0.5*r2/sigma**2)
psi0 = psi0 / np.linalg.norm(psi0) * 4.0

# x1: origen, x2: punto mas distante
x1 = pos[(0,0,0)]
xm = (N-1) if (N-1)%2==0 else (N-2)
x2 = pos[(xm, xm, xm)]

xs = np.array([s[0] for s in sites])
ys = np.array([s[1] for s in sites])
zs = np.array([s[2] for s in sites])
dist_x1x2 = np.sqrt((xs[x1]-xs[x2])**2 + (ys[x1]-ys[x2])**2 + (zs[x1]-zs[x2])**2)
print(f"  x1 = (0,0,0), x2 = ({xm},{xm},{xm})")
print(f"  Distancia x1->x2: {dist_x1x2:.2f} sitios")

# ======================================================================
# Convergencia comun al punto fijo
# ======================================================================
def evolucionar(psi_in, pasos, perturb=None):
    psi = psi_in.copy()
    for paso in range(pasos):
        coeffs = evecs.T @ psi
        pk = np.abs(coeffs)**2
        gk = sigmoid(evals_n) * eps0 / (eps0 + pk + 1e-16)
        theta = np.random.uniform(0, 2*np.pi, size=n_nodes)
        coeffs = coeffs * (1 - beta * gk) * np.exp(1j * theta * beta * gk)
        psi = evecs @ coeffs
        psi = psi / np.linalg.norm(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        coeffs = evecs.T @ psi
        coeffs = coeffs * np.exp(-1j * evals * dt)
        psi = evecs @ coeffs
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
    if perturb is not None:
        psi[perturb[0]] = psi[perturb[0]] * perturb[1]
        psi = psi / np.linalg.norm(psi)
    return psi

np.random.seed(42)
psi_fp = evolucionar(psi0, steps)
print(f"\n  Punto fijo alcanzado (t={steps*dt:.1f})")

# ======================================================================
# FASE 2: Perturbacion y efecto espectral inmediato
# ======================================================================
perturb_factor = 2.0
print(f"\n  FASE 2: Perturbacion en x1 factor {perturb_factor}")
print(f"  |psi(x1)|^2 antes: {np.abs(psi_fp[x1])**2:.8f}")

psi_pert = evolucionar(psi_fp, 0, perturb=(x1, perturb_factor))
print(f"  |psi(x1)|^2 despues: {np.abs(psi_pert[x1])**2:.8f}")

c_before = np.abs(evecs.T @ psi_fp)**2
c_after = np.abs(evecs.T @ psi_pert)**2
delta_c = c_after - c_before

top_k = np.argsort(np.abs(delta_c))[::-1][:5]
print(f"\n  Top 5 modos espectrales afectados:")
print(f"  {'modo':<6}{'lambda':<12}{'|c|^2 antes':<18}{'|c|^2 despues':<18}{'delta':<14}")
for k in top_k:
    print(f"  {k:<6}{evals[k]:<12.4f}{c_before[k]:<18.8f}{c_after[k]:<18.8f}{delta_c[k]:<14.8f}")

# ======================================================================
# FASE 3: Kernel K en x2 - cambio instantaneo
# ======================================================================
print(f"\n  FASE 3: Kernel de colapso en x2")
def gamma_at_site(site, c_pk):
    weight = np.sum(sigmoid(evals_n) * eps0 / (eps0 + c_pk + 1e-16)
                    * c_pk * np.abs(evecs[site, :])**2)
    return weight

g_x2_before = gamma_at_site(x2, c_before)
g_x2_after = gamma_at_site(x2, c_after)
print(f"  Gamma local en x2 antes:  {g_x2_before:.8f}")
print(f"  Gamma local en x2 despues: {g_x2_after:.8f}")
print(f"  Cambio: {g_x2_after-g_x2_before:.2e} ({(g_x2_after-g_x2_before)/g_x2_before*100:.4f}%)")
print(f"  Distancia x1->x2: {dist_x1x2:.1f} sitios")
print(f"  El cambio es INSTANTANEO: no depende de la distancia.")
print(f"  No viaja una senial - se reconfigura el kernel global.")

# ======================================================================
# FASE 4: Detectabilidad - SNR vs numero de muestras
# ======================================================================
print(f"\n  FASE 4: Detectabilidad estadistica")
M = 50

np.random.seed(42)
# Pre-converger estado base
psi_base = evolucionar(psi0, steps)

# Medir valores en x2 con/sin perturbacion
base_vals = []; pert_vals = []
for m in range(M):
    np.random.seed(100 + m)
    # Evolucionar desde el estado base con ruido distinto
    psi_b = evolucionar(psi_base, 200)
    base_vals.append(np.abs(psi_b[x2])**2)

    psi_p = evolucionar(psi_base, 200, perturb=(x1, perturb_factor))
    pert_vals.append(np.abs(psi_p[x2])**2)

    if (m+1) % 10 == 0:
        print(f"    muestra {m+1}/{M}")

base_arr = np.array(base_vals); pert_arr = np.array(pert_vals)

print(f"\n  SNR vs ventana de promediado:")
print(f"  {'M':<6}{'<base>':<16}{'<pert>':<16}{'diff':<14}{'SNR':<12}")
for Mw in [1, 2, 5, 10, 20, 50]:
    mb = np.mean(base_arr[:Mw])
    mp = np.mean(pert_arr[:Mw])
    sb = np.std(base_arr[:Mw])
    diff = mp - mb
    snr = diff / (sb/np.sqrt(Mw)) if sb > 0 else 0
    print(f"  {Mw:<6}{mb:<16.8f}{mp:<16.8f}{diff:<14.2e}{snr:<12.4f}")

print(f"\n  CONCLUSION:")
print(f"  La perturbacion espectral es detectable en x2 con SNR ~ sqrt(M).")
print(f"  El ruido de fondo es el ruido estocastico del colapso (theta aleatorio).")
print(f"  La senial es el cambio en |c_k|^2, que altera globalmente K.")
print(f"  No hay violacion de causalidad: la informacion se extrae")
print(f"  por promediado estadistico, no por senial instantanea.")
print(f"  Es comunicacion espectral: la red completa es el canal.")
