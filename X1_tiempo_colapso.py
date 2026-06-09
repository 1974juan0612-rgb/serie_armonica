import numpy as np, math

print("=" * 72)
print("  X1: COLAPSO EN K-SPACE  -  ESCALA TEMPORAL")
print("  Pregunta abierta: el colapso COMPLETA en t finito?")
print("=" * 72)

Z = 12; eps0 = 1.0/Z
eta = math.pi / (3*math.sqrt(2)); beta = (1 - eta) + eps0/2

# -----------------------------------------------------------------------
# Red FCC N=6, estado Bell en k0=0, k1=1
# Seguimiento de purificacion hasta t=20
# -----------------------------------------------------------------------
N = 6
sites = []; pos = {}; idx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (i + j + k) % 2 == 0:
                sites.append((i,j,k)); pos[(i,j,k)] = idx; idx += 1
n_nodes = len(sites)

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

psi_k = np.zeros(n_nodes, dtype=complex)
psi_k[0] = 1.0/math.sqrt(2); psi_k[1] = 1.0/math.sqrt(2)
psi = evecs @ psi_k
psi = psi / np.linalg.norm(psi) * 4.0

print(f"\n  Estado inicial: Bell en k0 (lambda=0) + k1 (lambda={evals[1]:.1f})")
print(f"  {'t':<8}{'ck0':<12}{'ck1':<12}{'S_espectral':<14}{'tasa_colapso':<16}{'t_estimado':<12}")
print(f"  " + "-" * 74)

dt = 0.005; n_pasos = 4000
t_anterior = 0; ck0_anterior = 0.5
t_estimado = float('inf')

for paso in range(0, n_pasos + 1, 50):
    t = paso * dt
    for _ in range(50):
        coeffs = evecs.T @ psi
        pk = np.abs(coeffs)**2
        gk = sigmoid(evals_n) * eps0 / (eps0 + pk + 1e-16)
        theta = np.random.uniform(0, 2*np.pi, size=n_nodes)
        coeffs = coeffs * (1 - beta * gk) * np.exp(1j * theta * beta * gk)
        psi = evecs @ coeffs; psi = psi / np.linalg.norm(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * Z * P * dt/2)
        coeffs = evecs.T @ psi
        coeffs = coeffs * np.exp(-1j * evals * dt)
        psi = evecs @ coeffs; psi = psi / np.linalg.norm(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * Z * P * dt/2)

    ck = np.abs(evecs.T @ psi)**2
    ck = ck / np.sum(ck)
    S = -np.sum(ck * np.log(ck + 1e-30))

    if t > 0 and ck[0] > ck0_anterior:
        tasa = (ck[0] - ck0_anterior) / (t - t_anterior)
        # Tiempo para alcanzar ck[0] = 1 (extrapolacion lineal)
        if tasa > 1e-8:
            t_estimado = t + (1.0 - ck[0]) / tasa
        else:
            t_estimado = float('inf')
    else:
        tasa = 0
        t_estimado = float('inf')

    ck0_anterior = ck[0]; t_anterior = t
    t_str = f"{t_estimado:.0f}" if t_estimado < 10000 else "inf"
    print(f"  {t:<8.2f}{ck[0]:<12.6f}{ck[1]:<12.6f}{S:<14.4f}{tasa:<16.8f}{t_str:<12}")

# -----------------------------------------------------------------------
# Ajuste: ck0(t) = 1 - A * exp(-t/tau)
# -----------------------------------------------------------------------
print(f"\n  AJUSTE EXPONENCIAL: ck0(t) = 1 - A * exp(-t/tau)")
print(f"  ===============================================")
# Recolectar datos
times = []; ck0s = []
psi_k = np.zeros(n_nodes, dtype=complex)
psi_k[0] = 1.0/math.sqrt(2); psi_k[1] = 1.0/math.sqrt(2)
psi = evecs @ psi_k; psi = psi / np.linalg.norm(psi) * 4.0

for paso in range(0, n_pasos + 1, 10):
    t = paso * dt
    for _ in range(10):
        coeffs = evecs.T @ psi
        pk = np.abs(coeffs)**2
        gk = sigmoid(evals_n) * eps0 / (eps0 + pk + 1e-16)
        theta = np.random.uniform(0, 2*np.pi, size=n_nodes)
        coeffs = coeffs * (1 - beta * gk) * np.exp(1j * theta * beta * gk)
        psi = evecs @ coeffs; psi = psi / np.linalg.norm(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * Z * P * dt/2)
        coeffs = evecs.T @ psi
        coeffs = coeffs * np.exp(-1j * evals * dt)
        psi = evecs @ coeffs; psi = psi / np.linalg.norm(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * Z * P * dt/2)
    ck = np.abs(evecs.T @ psi)**2; ck = ck / np.sum(ck)
    times.append(t); ck0s.append(ck[0])

# Ajuste exponencial a decaimiento 1-ck0(t)
ck0s = np.array(ck0s); times = np.array(times)
mask = times < 10  # primeros 10 unidades
y = 1.0 - ck0s[mask]
if np.any(y > 0):
    log_y = np.log(np.maximum(y, 1e-15))
    coeffs_poly = np.polyfit(times[mask], log_y, 1)
    tau = -1.0 / coeffs_poly[0]
    A = math.exp(coeffs_poly[1])
    print(f"  ck0(t) = 1 - {A:.4f} * exp(-t / {tau:.2f})")
    print(f"  Tiempo para ck0=0.95: {-tau * math.log((1-0.95)/A):.1f}")
    print(f"  Tiempo para ck0=0.99: {-tau * math.log((1-0.99)/A):.1f}")

# -----------------------------------------------------------------------
# Conclusion: pregunta abierta
# -----------------------------------------------------------------------
print(f"\n" + "=" * 72)
print("  PREGUNTA ABIERTA: colapso completo o asintotico?")
print("=" * 72)
print(f"""
  La tasa de colapso Gamma_k = sigmoid * eps0/(eps0+|ck|^2)
  se AUTOLIMITA: cuando |ck| es pequeno, Gamma_k es pequeno.

  Si la purificacion es asintotica (no completa en t finito):
  -> El universo SIEMPRE esta en transitorio
  -> Las galaxias JWST a z>10 son naturales
  -> Pero la CMB requiere homogeneidad a z~1100
  -> CONFLICTO si el transitorio es demasiado inhomogeneo

  Si la purificacion COMPLETA en t finito:
  -> El universo alcanza el punto fijo
  -> Se pierde el redshift espectral (z -> 0)
  -> Las galaxias JWST requieren otro mecanismo
  -> CONFLICTO con las observaciones de JWST

  Esta es la PREGUNTA CENTRAL que conecta A, B y C:
  la dinamica del colapso determina TODO el resto.

  Proximo paso: medir t_{1/2} para diferentes N y Normas,
  y determinar si tau es finito o infinito.
""")
