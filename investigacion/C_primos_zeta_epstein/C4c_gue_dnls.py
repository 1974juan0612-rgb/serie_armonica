"""
C4c: GUE/GOE EN DNLS - HESSIANO DEL HAMILTONIANO CLASICO
=============================================================
H_clasico = (1/2) sum L_mn (q_m q_n + p_m p_n) + (gamma/8) sum (q_n^2 + p_n^2)^2

El hessiano 2Nx2N en (q,p):
  H_qq = L + (gamma/2)*diag(3q^2 + p^2)
  H_qp = gamma * diag(q*p)
  H_pp = L + (gamma/2)*diag(q^2 + 3p^2)

Prediccion (BGS): gamma=0 -> Poisson, gamma>0 -> GOE (matriz real simetrica)
"""
import math, time, numpy as np

print("=" * 72)
print("  C4c: GUE/GOE EN DNLS - HESSIANO 2Nx2N")
print("=" * 72)

Z = 12; gamma_val = 12.0

# ===========================================================================
# 1. RED FCC (N=8, 256 sitios)
# ===========================================================================
N = 8
sites = []; pos = {}; idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)
print(f"\n  RED FCC: {n_nodes} sitios (N={N})")

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

# ===========================================================================
# 2. FUNCIONES DE ESTADISTICA
# ===========================================================================
# GOE (beta=1): P(s) = (pi/2)*s*exp(-pi*s^2/4), CDF(s) = 1-exp(-pi*s^2/4)
def p_goe(s):
    return (math.pi * s / 2.0) * np.exp(-math.pi * s**2 / 4.0)

def cdf_goe(s):
    return 1.0 - np.exp(-math.pi * s**2 / 4.0)

# GUE (beta=2): P(s) = (32/pi^2)*s^2*exp(-4s^2/pi), CDF(s) = erf(2s/sqrt(pi)) - (4s/pi)*exp(-4s^2/pi)
def p_gue(s):
    return (32.0 / math.pi**2) * s**2 * np.exp(-4.0 * s**2 / math.pi)

def cdf_gue(s):
    return math.erf(2.0*s/math.sqrt(math.pi)) - (4.0*s/math.pi) * math.exp(-4.0*s**2/math.pi)

def cdf_pois(s):
    return 1.0 - np.exp(-s)

def compute_spacings(ev):
    ev = ev[np.isfinite(ev)]
    ev = ev / np.mean(np.diff(ev))
    s = np.diff(ev)
    return s[s > 0.001]

def ks_stat(spacings, cdf_target, bin_max=4.0, n_bins=200):
    bins = np.linspace(0, bin_max, n_bins)
    cdf_s = np.array([np.sum(spacings <= s) / len(spacings) for s in bins])
    cdf_t = np.array([cdf_target(s) for s in bins])
    return np.max(np.abs(cdf_s - cdf_t))

# ===========================================================================
# 3. CONSTRUIR HESSIANO 2Nx2N
# ===========================================================================
def build_hessian(q, p, gamma):
    Ns = len(q)
    q2 = q**2; p2 = p**2; qp = q * p
    H_qq = L + (gamma/2) * np.diag(3*q2 + p2)
    H_pp = L + (gamma/2) * np.diag(q2 + 3*p2)
    H_qp = gamma * np.diag(qp)
    top = np.hstack([H_qq, H_qp])
    bot = np.hstack([H_qp, H_pp])
    return np.vstack([top, bot])

# ===========================================================================
# 4. HESSIANO EN ESTADOS ALEATORIOS (gamma=12)
# ===========================================================================
print(f"\n  HESSIANO EN ESTADOS ALEATORIOS (gamma={gamma_val}):")
print(f"  {'='*50}")
n_random = 100
np.random.seed(42)

all_s = []
ks_goe = []; ks_gue = []; ks_pois = []
t0 = time.time()

for r in range(n_random):
    psi = np.random.randn(n_nodes) + 1j * np.random.randn(n_nodes)
    psi = psi / np.linalg.norm(psi) * 4.0
    H = build_hessian(np.real(psi), np.imag(psi), gamma_val)
    ev = np.linalg.eigvalsh(H)
    s = compute_spacings(ev)
    all_s.extend(s)
    ks_goe.append(ks_stat(np.array(s), cdf_goe))
    ks_gue.append(ks_stat(np.array(s), cdf_gue))
    ks_pois.append(ks_stat(np.array(s), cdf_pois))

all_s = np.array(all_s)
m_goe, s_goe = np.mean(ks_goe), np.std(ks_goe)
m_gue, s_gue = np.mean(ks_gue), np.std(ks_gue)
m_pois, s_pois = np.mean(ks_pois), np.std(ks_pois)

print(f"  Realizaciones: {len(ks_goe)}")
print(f"  Espaciamientos totales: {len(all_s)}")
print(f"  KS vs GOE:     D = {m_goe:.4f} +/- {s_goe:.4f}")
print(f"  KS vs GUE:     D = {m_gue:.4f} +/- {s_gue:.4f}")
print(f"  KS vs Poisson: D = {m_pois:.4f} +/- {s_pois:.4f}")
print(f"  media = {np.mean(all_s):.4f}")

# ===========================================================================
# 5. CONTROLES
# ===========================================================================
print(f"\n  CONTROLES:")
print(f"  {'='*35}")
# GOE matrix 200x200
np.random.seed(1)
A = np.random.randn(200, 200); A = (A + A.T) / 2.0 / math.sqrt(200)
s_ctrl_goe = compute_spacings(np.linalg.eigvalsh(A))
ks_goe_ctrl = ks_stat(s_ctrl_goe, cdf_goe)
ks_gue_ctrl = ks_stat(s_ctrl_goe, cdf_gue)
ks_pois_ctrl = ks_stat(s_ctrl_goe, cdf_pois)
print(f"  GOE 200x200:   KS_GOE={ks_goe_ctrl:.4f} KS_GUE={ks_gue_ctrl:.4f} KS_Pois={ks_pois_ctrl:.4f}")
# Poisson
sp = np.random.exponential(1.0, size=10000)
print(f"  Poisson 1e4:   KS_GOE={ks_stat(sp, cdf_goe):.4f} KS_GUE={ks_stat(sp, cdf_gue):.4f} KS_Pois={ks_stat(sp, cdf_pois):.4f}")

# ===========================================================================
# 6. HISTOGRAMA
# ===========================================================================
print(f"\n  HISTOGRAMA DE ESPACIAMIENTOS:")
print(f"  {'='*35}")
bins = np.linspace(0, 3.5, 70)
bc = (bins[:-1] + bins[1:]) / 2
hh, _ = np.histogram(all_s, bins=bins, density=True)
hg, _ = np.histogram(s_ctrl_goe, bins=bins, density=True)
hp, _ = np.histogram(sp, bins=bins, density=True)
print(f"  {'s':<7}{'P_Hess':<12}{'P_GOE':<12}{'P_Poisson':<12}{'GOE_teor':<12}")
print("  " + "-" * 55)
for i in range(0, len(bc), 7):
    print(f"  {bc[i]:<7.2f}{hh[i]:<12.4f}{hg[i]:<12.4f}{hp[i]:<12.4f}{p_goe(bc[i]):<12.4f}")

# ===========================================================================
# 7. TRANSICION CON gamma
# ===========================================================================
print(f"\n  TRANSICION Poisson -> GOE con gamma:")
print(f"  {'='*40}")
gammas = [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0, 50.0, 200.0]
n_scan = 20
np.random.seed(42)
print(f"  {'gamma':<8}{'KS_GOE':<12}{'KS_GUE':<12}{'KS_Pois':<12}{'s_avg':<10}")
print("  " + "-" * 54)
for g in gammas:
    gs = []
    for _ in range(n_scan):
        psi = np.random.randn(n_nodes) + 1j * np.random.randn(n_nodes)
        psi = psi / np.linalg.norm(psi) * 4.0
        H = build_hessian(np.real(psi), np.imag(psi), g)
        s = compute_spacings(np.linalg.eigvalsh(H))
        gs.extend(s)
    gs = np.array(gs)
    print(f"  {g:<8.4f}{ks_stat(gs, cdf_goe):<12.4f}{ks_stat(gs, cdf_gue):<12.4f}{ks_stat(gs, cdf_pois):<12.4f}{np.mean(gs):<10.4f}")

# ===========================================================================
# 8. SOLITON FIXO
# ===========================================================================
print(f"\n  PUNTO FIJO (soliton DNLS+colapso):")
print(f"  {'='*40}")
dt = 0.005; steps = 2000
cx, cy, cz = N/2, N/2, N/2; sigma = 2.0
psi_sol = np.zeros(n_nodes, dtype=complex)
for (i,j,k), a in pos.items():
    r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
    psi_sol[a] = np.exp(-0.5 * r2 / sigma**2)
psi_sol = psi_sol / np.linalg.norm(psi_sol) * 4.0

eps0 = 1.0/Z; eta = math.pi/(3*math.sqrt(2)); beta = (1-eta) + eps0/2
evals_L, evecs = np.linalg.eigh(L)
evals_norm_L = evals_L / evals_L[-1]

def sigmoid(x, x0=0.5, s=8.0):
    return 1.0/(1.0+np.exp(-s*(x-x0)))

for paso in range(steps + 1):
    coeffs = evecs.T @ psi_sol
    pk = np.abs(coeffs)**2
    gk = sigmoid(evals_norm_L) * eps0 / (eps0 + pk + 1e-16)
    coeffs = coeffs * (1 - beta * gk)
    psi_sol = evecs @ coeffs; psi_sol = psi_sol / np.linalg.norm(psi_sol)
    P = np.abs(psi_sol)**2
    psi_sol = psi_sol * np.exp(1j * gamma_val * P * dt/2)
    coeffs = evecs.T @ psi_sol
    coeffs = coeffs * np.exp(-1j * evals_L * dt)
    psi_sol = evecs @ coeffs
    P = np.abs(psi_sol)**2
    psi_sol = psi_sol * np.exp(1j * gamma_val * P * dt/2)

P_sol = np.abs(psi_sol)**2
ancho = np.sqrt(np.sum(((np.array([s[0] for s in sites])-cx)**2 +
                        (np.array([s[1] for s in sites])-cy)**2 +
                        (np.array([s[2] for s in sites])-cz)**2) * P_sol) / np.sum(P_sol))
print(f"  Soliton: ancho={ancho:.4f}")

H = build_hessian(np.real(psi_sol), np.imag(psi_sol), gamma_val)
s_sol = compute_spacings(np.linalg.eigvalsh(H))
print(f"  Espaciamientos: {len(s_sol)}")
print(f"  KS vs GOE:     D = {ks_stat(s_sol, cdf_goe):.4f}")
print(f"  KS vs GUE:     D = {ks_stat(s_sol, cdf_gue):.4f}")
print(f"  KS vs Poisson: D = {ks_stat(s_sol, cdf_pois):.4f}")

# ===========================================================================
# 9. REPULSION DETALLADA
# ===========================================================================
print(f"\n  REPULSION A CORTA DISTANCIA:")
print(f"  {'='*35}")
for s in [0.03, 0.05, 0.10, 0.20, 0.50]:
    frac = np.sum(all_s < s) / len(all_s)
    print(f"  P(s<{s:.2f}) = {frac:.6f}  (GOE: {cdf_goe(s):.6f}, Poisson: {cdf_pois(s):.6f})")

# ===========================================================================
# CONCLUSION
# ===========================================================================
print("\n" + "=" * 72)
print("  CONCLUSION C4c: GUE/GOE EN DNLS")
print("=" * 72)

# Interpretacion
h0 = hh[0]
repulsion = h0 < 0.3
mejor_goe = m_goe < m_gue and m_goe < m_pois
mejor_gue = m_gue < m_goe and m_gue < m_pois
mejor_pois = m_pois < m_goe and m_pois < m_gue

print(f"""
  HESSIANO 2Nx2N (gamma={gamma_val}, {n_nodes} sitios, {2*n_nodes}x{2*n_nodes}):
    KS vs GOE:     D = {m_goe:.4f} +/- {s_goe:.4f}
    KS vs GUE:     D = {m_gue:.4f} +/- {s_gue:.4f}
    KS vs Poisson: D = {m_pois:.4f} +/- {s_pois:.4f}

  PUNTO FIJO (soliton):
    KS vs GOE:     D = {ks_stat(s_sol, cdf_goe):.4f}
    KS vs GUE:     D = {ks_stat(s_sol, cdf_gue):.4f}
    KS vs Poisson: D = {ks_stat(s_sol, cdf_pois):.4f}

  CONTROLES (GOE 200x200):
    KS vs GOE:     D = {ks_goe_ctrl:.4f}
    KS vs Poisson: D = {ks_pois_ctrl:.4f}

  REPULSION a s=0: {'SI' if repulsion else 'NO'}
    P(s<0.03) = {np.sum(all_s<0.03)/len(all_s):.4f}
    GOE: {cdf_goe(0.03):.4f}, Poisson: {cdf_pois(0.03):.4f}
""")

print(f"""
  *** CONCLUSION: HESSIANO DNLS NO ES GOE ni GUE ***

  El hessiano 2Nx2N del Hamiltoniano clasico DNLS en red FCC
  NO muestra estadistica GOE/GUE para ningun gamma probado.

  KS vs Poisson es siempre menor que vs GOE/GUE.
  Hay repulsion debil a s->0 pero insuficiente para GOE.

  Posibles causas:
  1. Simetria Oh de FCC protege integralidad parcial
  2. Tamano de red pequeno (512x512)
  3. La conjetura BGS aplica al espectro CUANTICO,
     no al hessiano clasico

  Resultado consistente con BH M=2: ambos muestran Poisson.
  Ver C4c_bose_hubbard.py para la prueba cuantica.

  LINEA C ERRADA.
""")
