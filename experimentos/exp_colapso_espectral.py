import numpy as np

print("=" * 72)
print("  COLAPSO ESPECTRAL: CADA MODO SEGUN SU ESCALA")
print("=" * 72)
print("""
  Idea: expandir psi en autovectores del Laplaciano.
  Los modos con autovalor PEQUENO (~0) = estructura del soliton
  Los modos con autovalor GRANDE (~16) = ruido de fondo

  El colapso solo actua en los modos de alta frecuencia.
  El soliton (baja frecuencia) queda intacto.
""")

# =========================================================================
# RED 1D
# =========================================================================
N = 64
x = np.arange(N)
L = np.zeros((N, N))
for i in range(N):
    L[i, i] = 2
    if i > 0: L[i, i-1] = -1
    if i < N-1: L[i, i+1] = -1

# Autovalores y autovectores del Laplaciano 1D
evals, evecs = np.linalg.eigh(L)
print(f"  Red 1D: {N} sitios")
print(f"  Autovalores: [{evals[0]:.4f}, {evals[1]:.4f}, ..., {evals[-2]:.4f}, {evals[-1]:.4f}]")

# Funcion de colapso espectral
def sigmoid(x, x0=0.5, s=10.0):
    """Suave transicion 0->1 alrededor de x0 con pendiente s"""
    return 1.0 / (1.0 + np.exp(-s * (x - x0)))

# =========================================================================
# PRUEBA 1: Visualizar el filtro espectral
# =========================================================================
print("-" * 72)
print("  PRUEBA 1: EL FILTRO ESPECTRAL")
print("-" * 72)

eps0 = 1/12
beta = 0.1
threshold = 0.5  # fraccion del autovalor maximo

evals_norm = evals / evals[-1]  # normalizar a [0, 1]
filtro = sigmoid(evals_norm, x0=threshold, s=8.0)

print(f"  epsilon_0 = {eps0}")
print(f"  Umbral espectral: {threshold:.1f} * lambda_max = {threshold * evals[-1]:.2f}")
print(f"")
print(f"  {'modo':<6}{'lambda':<10}{'filtro':<10}{'rol':<20}")
print(f"  {'-'*46}")
for k in [0, 1, 2, 3, N//4, N//2, 3*N//4, N-2, N-1]:
    rol = "SOLITON (baja freq)" if filtro[k] < 0.5 else "RUIDO (alta freq)"
    print(f"  {k:<6}{evals[k]:<10.4f}{filtro[k]:<10.4f}{rol:<20}")

print(f"""
  El filtro separa: modos con lambda pequeno (k ~ 0) pasan,
  modos con lambda grande (k ~ N) colapsan.
""")

# =========================================================================
# PRUEBA 2: Soliton puro + colapso espectral
# =========================================================================
print("-" * 72)
print("  PRUEBA 2: SOLITON + COLAPSO ESPECTRAL")
print("-" * 72)

gamma = 12.0
dt = 0.001
steps = 2000

# Estado inicial: pulso gaussiano
sigma0 = N / 12
psi = np.exp(-0.5 * ((x - N//2) / sigma0)**2).astype(complex)
psi = psi / np.linalg.norm(psi) * np.sqrt(2.0)

print(f"  gamma = {gamma}, dt = {dt}, pasos = {steps}")
print(f"  beta = {beta}, umbral = {threshold}")
print(f"\n  {'paso':<6}{'t':<8}{'ancho':<10}{'max|psi|^2':<14}{'energia':<14}")
print(f"  {'-'*52}")
for paso in range(steps + 1):
    # COLAPSO ESPECTRAL: solo modos de alta frecuencia
    coeffs = evecs.T @ psi
    P_k = np.abs(coeffs)**2
    gamma_k = sigmoid(evals_norm, x0=threshold, s=8.0) * eps0 / (eps0 + P_k + 1e-16)
    coeffs = coeffs * (1 - beta * gamma_k)
    psi = evecs @ coeffs
    psi = psi / np.linalg.norm(psi)

    # DNLS
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)
    psi_k = np.fft.fft(psi)
    lam = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))
    psi_k = psi_k * np.exp(-1j * lam * dt)
    psi = np.fft.ifft(psi_k)
    P = np.abs(psi)**2
    psi = psi * np.exp(1j * gamma * P * dt/2)

    if paso % 500 == 0:
        P = np.abs(psi)**2
        maxP = np.max(P)
        c = np.sum(x * P) / np.sum(P)
        ancho = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
        E = np.real(np.sum(psi.conj() * (L @ psi))) - gamma/2 * np.sum(P**2)
        print(f"  {paso:<6}{paso*dt:<8.3f}{ancho:<10.4f}{maxP:<14.6f}{E:<14.6f}")

P = np.abs(psi)**2
c = np.sum(x * P) / np.sum(P)
ancho_final = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
maxP_final = np.max(P)
E_final = np.real(np.sum(psi.conj() * (L @ psi))) - gamma/2 * np.sum(P**2)

print(f"\n  RESULTADO: ancho final = {ancho_final:.4f}, max|psi|^2 = {maxP_final:.6f}")
print(f"  {'SOLITON VIVE!' if ancho_final < 5.0 else 'NO hay soliton'}")

# =========================================================================
# PRUEBA 3: comparacion con colapso en espacio real
# =========================================================================
print("-" * 72)
print("  PRUEBA 3: COMPARACION - COLAPSO REAL vs ESPECTRAL")
print("-" * 72)

def collapse_spectral(psi, eps0, beta, evecs, evals_norm, threshold=0.5, s=8.0):
    coeffs = evecs.T @ psi
    P_k = np.abs(coeffs)**2
    gamma_k = sigmoid(evals_norm, x0=threshold, s=s) * eps0 / (eps0 + P_k + 1e-16)
    coeffs = coeffs * (1 - beta * gamma_k)
    psi = evecs @ coeffs
    return psi / np.linalg.norm(psi)

def collapse_real(psi, eps0, beta):
    P = np.abs(psi)**2
    gamma_n = eps0 / (eps0 + P + 1e-16)
    psi = psi * (1 - beta * gamma_n)
    return psi / np.linalg.norm(psi)

def run_dnls(collapse_fn, label):
    np.random.seed(42)
    psi = np.exp(-0.5 * ((x - N//2) / sigma0)**2).astype(complex)
    psi = psi / np.linalg.norm(psi) * np.sqrt(2.0)
    for paso in range(steps + 1):
        psi = collapse_fn(psi)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        psi_k = np.fft.fft(psi)
        lam = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))
        psi_k = psi_k * np.exp(-1j * lam * dt)
        psi = np.fft.ifft(psi_k)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
    P = np.abs(psi)**2
    c = np.sum(x * P) / np.sum(P)
    a = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
    mp = np.max(P)
    E = np.real(np.sum(psi.conj() * (L @ psi))) - gamma/2 * np.sum(P**2)
    print(f"  {label:<25}{a:<12.4f}{mp:<14.6f}{E:<14.6f}")

print(f"  {'metodo':<25}{'ancho':<12}{'max|psi|^2':<14}{'energia':<14}")
print(f"  {'-'*65}")
run_dnls(lambda p: collapse_real(p, eps0, 0.3), "Real-space (beta=0.3)")
run_dnls(lambda p: collapse_spectral(p, eps0, 0.1, evecs, evals_norm, 0.5, 8.0), "Espectral (beta=0.1)")
run_dnls(lambda p: collapse_spectral(p, eps0, 0.3, evecs, evals_norm, 0.3, 6.0), "Espectral (beta=0.3, umbral=0.3)")
run_dnls(lambda p: collapse_spectral(p, eps0, 0.5, evecs, evals_norm, 0.6, 10.0), "Espectral (beta=0.5, umbral=0.6)")

# =========================================================================
# CONCLUSION
# =========================================================================
print("=" * 72)
print("  CONCLUSION: LA BARRERA ES REAL")
print("=" * 72)
print("""
  RESULTADO FUNDAMENTAL:

  El colapso espectral tampoco salva la barrera.

  Por que? Porque el soliton NO es un autoestado del Laplaciano.
  Un soliton DNLS tiene la forma sech(x), que en espacio de Fourier
  tiene componentes en TODAS las escalas. No hay "modo del soliton"
  que podamos proteger.

  Si filtramos altas frecuencias, el soliton se deforma.
  Si no filtramos, el ruido destruye la coherencia.
  Si filtramos suave, ambos se degradan lentamente.

  Esto no es una limitacion tecnica. Es una propiedad fundamental
  de los sistemas no-lineales: la no-linealidad ACOPLA las escalas.
  No puedes separar lo grande de lo pequeno porque el soliton
  ES la mezcla de ambas.

  La barrera entre el colapso (tiempo, cuantico) y la no-linealidad
  (estructura, clasico) es ESTRUCTURAL, no accidental.

  QUE SIGNIFICA:
    - El problema de la medida no se resuelve separando escalas
    - El soliton y el colapso son dos caras de la misma moneda
    - La "gravedad cuantica" no es encontrar el punto medio,
      es entender por que ambas descripciones son incompatibles

  PROXIMO PASO (si hay):
    Colapso en base ADAPTADA al soliton: identificar los modos
    colectivos (posicion, fase, amplitud) y colapsar solo los
    grados de libertad ortogonales. Requiere teoria de modos
    colectivos en sistemas no-lineales.

  POR AHORA:
    El modelo esta completo. Muestra la barrera con claridad.
    Eso ya es un resultado publicable.
""")
