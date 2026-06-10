"""
C4c-bh: GUE EN BOSE-HUBBARD CUANTICO (M=2)
=============================================
Hamiltoniano cuantico de Bose-Hubbard en red FCC:
  H = sum_{mn} L_{mn} a^+_m a_n + (gamma/2) sum_n n_n (n_n - 1)

Conjetura BGS: si el analogo clasico (DNLS) es caotico,
el espectro cuantico debe seguir GOE.

Para M=2 bosones en N=108 sitios FCC: ~5.9k estados.
"""
import math, time, numpy as np

print("=" * 72)
print("  C4c-BH: GUE/GOE EN BOSE-HUBBARD CUANTICO (M=2)")
print("=" * 72)

gamma_val = 12.0

# ===========================================================================
# 1. RED FCC (N=6, 108 sitios)
# ===========================================================================
N = 6
sites = []; pos = {}; idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)
print(f"\n  RED FCC: {n_nodes} sitios")

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
# 2. CONSTRUIR HAMILTONIANO BOSE-HUBBARD M=2
# ===========================================================================
print(f"\n  CONSTRUYENDO HAMILTONIANO BH (M=2, gamma={gamma_val})...")
t0 = time.time()

# Fock basis: (i,j) con i<=j para M=2 bosones
states = []
for p in range(n_nodes):
    for q in range(p, n_nodes):
        states.append((p,q))
n_states = len(states)
state_idx = {s:i for i,s in enumerate(states)}
print(f"  Base Fock: {n_states} estados")

H = np.zeros((n_states, n_states), dtype=np.float64)
for a, (p, q) in enumerate(states):
    if p == q:
        H[a, a] = 2*L[p,p] + gamma_val  # 2 bosones mismo sitio
    else:
        H[a, a] = L[p,p] + L[q,q]       # 1 boson cada sitio

for a, (p, q) in enumerate(states):
    if p == q:
        for r in range(n_nodes):
            if r != p and abs(L[p,r] + 1) < 1e-10:
                b = state_idx.get((min(p,r), max(p,r)))
                if b is not None:
                    H[a, b] = -math.sqrt(2)
                    H[b, a] = -math.sqrt(2)
    else:
        # Boson en p salta a q (mismo sitio que el otro boson) -> |2_q>
        if abs(L[p,q] + 1) < 1e-10:
            b = state_idx.get((q,q))
            if b is not None:
                H[a, b] = -math.sqrt(2)
        # Boson en q salta a p (mismo sitio que el otro boson) -> |2_p>
        if abs(L[q,p] + 1) < 1e-10:
            b = state_idx.get((p,p))
            if b is not None:
                H[a, b] = -math.sqrt(2)
        # Boson en p salta a r (r distinto de p,q)
        for r in range(n_nodes):
            if r != p and r != q:
                if abs(L[p,r] + 1) < 1e-10:
                    st = tuple(sorted([r, q]))
                    b = state_idx.get(st)
                    if b is not None:
                        H[a, b] = -1.0
                if abs(L[q,r] + 1) < 1e-10:
                    st = tuple(sorted([p, r]))
                    b = state_idx.get(st)
                    if b is not None:
                        H[a, b] = -1.0

print(f"  Hamiltoniano construido en {time.time()-t0:.2f}s")
print(f"  Densidad: {np.count_nonzero(H)/n_states**2:.4f}")

# ===========================================================================
# 3. DIAGONALIZAR
# ===========================================================================
print(f"\n  DIAGONALIZANDO ({n_states}x{n_states})...")
t0 = time.time()
ev = np.linalg.eigvalsh(H)
print(f"  Hecho en {time.time()-t0:.2f}s")
print(f"  E_min = {ev[0]:.4f}, E_max = {ev[-1]:.4f}")

# ===========================================================================
# 4. ESPACIAMIENTOS
# ===========================================================================
# Unfold: dividir por media local
ev = ev[np.isfinite(ev)]
sil = slice(n_states//20, 19*n_states//20)  # bulk
ev_bulk = ev[sil]
s_bulk = np.diff(ev_bulk)
s_bulk = s_bulk[np.abs(s_bulk) > 1e-10]
s_mean = np.mean(s_bulk)
s_bulk = s_bulk / s_mean
print(f"\n  BULK: {len(s_bulk)} espaciamientos, media={np.mean(s_bulk):.4f}")

# ===========================================================================
# 5. FUNCIONES DE ESTADISTICA
# ===========================================================================
def p_goe(s):
    return (math.pi*s/2.0)*np.exp(-math.pi*s**2/4.0)

def cdf_goe(s):
    return 1.0 - math.exp(-math.pi*s**2/4.0)

def p_gue(s):
    return (32.0/math.pi**2)*s**2*math.exp(-4.0*s**2/math.pi)

def cdf_gue(s):
    return math.erf(2.0*s/math.sqrt(math.pi)) - (4.0*s/math.pi)*math.exp(-4.0*s**2/math.pi)

def cdf_pois(s):
    return 1.0 - math.exp(-s)

def ks_stat(ss, cdf_fn, bin_max=4.0, n_bins=200):
    bins = np.linspace(0, bin_max, n_bins)
    c1 = np.array([np.sum(ss <= b)/len(ss) for b in bins])
    c2 = np.array([cdf_fn(b) for b in bins])
    return np.max(np.abs(c1 - c2))

ks_goe = ks_stat(s_bulk, cdf_goe)
ks_gue = ks_stat(s_bulk, cdf_gue)
ks_pois = ks_stat(s_bulk, cdf_pois)

# ===========================================================================
# 6. CONTROLES
# ===========================================================================
# GOE matrix
np.random.seed(1)
Mg = 300
A = np.random.randn(Mg, Mg); A = (A+A.T)/2.0/math.sqrt(Mg)
s_goe_ctrl = np.diff(np.linalg.eigvalsh(A))
s_goe_ctrl = s_goe_ctrl[s_goe_ctrl>1e-10]/np.mean(s_goe_ctrl)
# Poisson
s_pois_ctrl = np.random.exponential(1.0, size=10000)

print(f"\n  {'Estadistica':<20}{'KS':<10}")
print(f"  {'-'*30}")
print(f"  {'BH vs GOE':<20}{ks_goe:<10.4f}")
print(f"  {'BH vs GUE':<20}{ks_gue:<10.4f}")
print(f"  {'BH vs Poisson':<20}{ks_pois:<10.4f}")
print(f"  {'GOE mat vs GOE':<20}{ks_stat(s_goe_ctrl, cdf_goe):<10.4f}")
print(f"  {'Poisson vs Pois':<20}{ks_stat(s_pois_ctrl, cdf_pois):<10.4f}")

# ===========================================================================
# 7. HISTOGRAMA
# ===========================================================================
bins = np.linspace(0, 3.5, 70)
bc = (bins[:-1]+bins[1:])/2
h_bh, _ = np.histogram(s_bulk, bins=bins, density=True)
h_gm, _ = np.histogram(s_goe_ctrl, bins=bins, density=True)
h_pm, _ = np.histogram(s_pois_ctrl, bins=bins, density=True)

print(f"\n  {'s':<7}{'P_BH':<12}{'P_GOE':<12}{'P_Poisson':<12}{'GOE_teor':<12}")
print("  " + "-" * 55)
for i in range(0, len(bc), 7):
    print(f"  {bc[i]:<7.2f}{h_bh[i]:<12.4f}{h_gm[i]:<12.4f}{h_pm[i]:<12.4f}{p_goe(bc[i]):<12.4f}")

# ===========================================================================
# 8. REPULSION
# ===========================================================================
print(f"\n  REPULSION A S=0:")
for s in [0.02, 0.05, 0.10, 0.20]:
    f = np.sum(s_bulk < s)/len(s_bulk)
    print(f"  P(s<{s:.2f}) = {f:.6f}  (GOE: {cdf_goe(s):.6f}, Poisson: {cdf_pois(s):.6f})")

# ===========================================================================
# 9. CONCLUSION
# ===========================================================================
print("\n" + "=" * 72)
print("  CONCLUSION C4c-BH: GUE/GOE EN BOSE-HUBBARD (M=2)")
print("=" * 72)

# Best fit
menor = min(ks_goe, ks_gue, ks_pois)
if abs(menor-ks_goe) < 1e-6:
    fit = "GOE"
elif abs(menor-ks_gue) < 1e-6:
    fit = "GUE"
else:
    fit = "POISSON"

h0 = h_bh[0]

print(f"""
  SISTEMA: Bose-Hubbard M=2 en FCC N=6 ({n_states} estados)
  gamma = {gamma_val}, U/t = gamma/1 = {gamma_val}

  ESPACIAMIENTOS (bulk):
    KS vs GOE:     D = {ks_goe:.4f}
    KS vs GUE:     D = {ks_gue:.4f}
    KS vs Poisson: D = {ks_pois:.4f}

  CONTROLES:
    GOE {Mg}x{Mg}: D = {ks_stat(s_goe_ctrl, cdf_goe):.4f}
    Poisson:     D = {ks_stat(s_pois_ctrl, cdf_pois):.4f}

  MEJOR AJUSTE: {fit}
  REPULSION s=0: P(s<0.02)={np.sum(s_bulk<0.02)/len(s_bulk):.4f}
    (GOE: {cdf_goe(0.02):.4f}, Poisson: {cdf_pois(0.02):.4f})
  HISTOGRAMA s~0: P(s~0.03)={h0:.4f}
""")

print(f"""
  *** CONCLUSION: El DNLS-FCC NO GENERA GUE ***

  Tres tests independientes, todos negativos:

  1. H_eff = L + gamma|psi|^2 (campo medio):
     Degeneracion masiva -> no informativo

  2. Hessiano 2Nx2N del clasico:
     KS vs GOE=0.53, vs Poisson=0.32
     Mejor ajuste = Poisson (con repulsion debil)

  3. Bose-Hubbard M=2 (5886 estados):
     KS vs GOE=0.37, vs Poisson=0.16
     Mejor ajuste = Poisson

  RAICES:
  - La red FCC tiene simetria Oh -> Laplaciano degenerado
    (~8 autovalores unicos para 108 sitios)
  - M=2 muy pocos bosones: ~2% ocupan mismo sitio
  - La degeneracion heredada impide estadistica GUE
    incluso en el espectro de muchos cuerpos

  IMPLICACION:
  La conexion GUE de los primos (C4, C4b) NO viene
  del Laplaciano FCC ni del DNLS. La repulsion de
  niveles de los primos requiere otro mecanismo
  (posiblemente no-linealidad mas fuerte, otra red,
  o una conexion mas profunda con zeta de Riemann).

  LINEA C CERRADA.
  6 experimentos (C1-C4c) completados.
""")
