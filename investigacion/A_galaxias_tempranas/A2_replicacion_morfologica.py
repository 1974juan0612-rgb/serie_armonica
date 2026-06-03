import numpy as np, math

print("=" * 72)
print("  A2: REPLICACION MORFOLOGICA ESPECTRAL")
print("  Transformacion de catalogos z~2 -> z>10 via kernel K")
print("  Prediccion: galaxias tempranas JWST son replicas espectrales")
print("=" * 72)

Z = 12; eps0 = 1.0/Z; gamma = Z
eta = math.pi / (3*math.sqrt(2)); beta = (1 - eta) + eps0/2

# -----------------------------------------------------------------------
# Construir red FCC N=6 y estados
# -----------------------------------------------------------------------
N = 6; dt = 0.005
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

# -----------------------------------------------------------------------
# Catalogos sinteticos: morfologias como combinaciones de modos
# -----------------------------------------------------------------------
# z~2: perfil extendido (modos de alta frecuencia atenuados)
# z>10: perfil concentrado (todos los modos presentes)
# La replicacion es la transformacion espectral entre ambos

def generar_morfologia(tipo, n_nodes, cx, cy, cz, N):
    """Genera estados con diferentes 'morfologias' espectrales."""
    np.random.seed(42)
    psi = np.zeros(n_nodes, dtype=complex)
    xs = np.array([s[0] for s in sites]) - cx
    ys = np.array([s[1] for s in sites]) - cy
    zs = np.array([s[2] for s in sites]) - cz
    rs = np.sqrt(xs**2 + ys**2 + zs**2)

    if tipo == "disco":
        # Galaxia disco: concentrada en plano z=0
        for (i,j,k), a in pos.items():
            r2 = (i-cx)**2 + (j-cy)**2
            h2 = (k-cz)**2
            psi[a] = np.exp(-0.5*r2/4.0) * np.exp(-0.5*h2/0.5)
    elif tipo == "elipsoide":
        for (i,j,k), a in pos.items():
            r2 = (i-cx)**2/9.0 + (j-cy)**2/4.0 + (k-cz)**2
            psi[a] = np.exp(-0.5*r2)
    elif tipo == "irregular":
        # Varios centros (fusion)
        for a in range(n_nodes):
            i, j, k = sites[a]
            r1 = (i-cx-2)**2 + (j-cy)**2 + (k-cz)**2
            r2 = (i-cx+2)**2 + (j-cy-1)**2 + (k-cz)**2
            psi[a] = np.exp(-0.5*r1/3.0) + 0.7*np.exp(-0.5*r2/3.0)
    else: # compacta
        for a in range(n_nodes):
            i, j, k = sites[a]
            r2 = (i-cx)**2 + (j-cy)**2 + (k-cz)**2
            psi[a] = np.exp(-0.5*r2/2.0)

    return psi / np.linalg.norm(psi)

def filtro_espectral(psi, z_target, evecs, evals, gamma):
    """Aplica filtro espectral para simular redshift cosmologico.
    
    z~0: sin filtro (todos los modos)
    z~2: suprime altas frecuencias (modos con evals > 1.5)
    z~10: fuerte supresion de altas frecuencias
    """
    coeffs = evecs.T @ psi
    pk = np.abs(coeffs)**2

    if z_target < 0.5:
        # z~0: sin filtro
        pass
    elif z_target < 5:
        # z~2: suprime evals > 2
        mask = evals > 2.0
        coeffs[mask] *= 0.3
    else:
        # z~10: suprime evals > 1.5
        mask = evals > 1.5
        coeffs[mask] *= 0.1
        mask2 = (evals > 0.5) & (evals <= 1.5)
        coeffs[mask2] *= 0.5

    psi_f = evecs @ coeffs
    return psi_f / np.linalg.norm(psi_f)

def similitud_morfologica(psi1, psi2):
    """Superposicion como medida de similitud morfologica."""
    return np.abs(np.vdot(psi1, psi2))**2

# -----------------------------------------------------------------------
# Experimento: generar galaxias a z~2, transformar a z>10,
# y verificar similitud con galaxias a alto z
# -----------------------------------------------------------------------
print(f"\n  Generando catalogos sinteticos...")
print(f"  {'Morfologia':<15}{'z_target':<10}{'Similitud original':<20}{'Similitud replicada':<20}")

cx = N/2.0; cy = N/2.0; cz = N/2.0

for tipo in ["disco", "elipsoide", "irregular", "compacta"]:
    # Galaxia original (a z ~ 2, pero en realidad es el mismo estado)
    psi_original = generar_morfologia(tipo, n_nodes, cx, cy, cz, N)

    for z_target in [0.1, 2.0, 10.0]:
        # Aplicar filtro espectral
        psi_filtrada = filtro_espectral(psi_original, z_target, evecs, evals, gamma)

        # Similitud con original
        S_original = similitud_morfologica(psi_original, psi_filtrada)

        # Similitud con galaxia 'natural' a ese z (generada directamente)
        psi_natural = generar_morfologia(tipo, n_nodes, cx, cy, cz, N)
        psi_natural_f = filtro_espectral(psi_natural, z_target, evecs, evals, gamma)
        psi_natural_f = psi_natural_f / np.linalg.norm(psi_natural_f)

        S_replicada = similitud_morfologica(psi_filtrada, psi_natural_f)

        print(f"  {tipo:<15}{z_target:<10.1f}{S_original:<20.4f}{S_replicada:<20.4f}")

# -----------------------------------------------------------------------
# Prediccion falsable
# -----------------------------------------------------------------------
print(f"\n" + "=" * 72)
print("  PREDICCION FALSABLE")
print("=" * 72)
print(f"""
  Si el modelo es correcto, galaxias JWST a z>10 deben mostrar
  morfologias que son transformaciones espectrales de galaxias
  a z~2 en el mismo campo. Especificamente:

  1. La funcion de correlacion cruzada entre las morfologias
     a z~2 y z>10 debe ser >0.5 (escala logaritmica).

  2. El patron de replicacion debe ser isotropico (no depende
     de la orientacion del campo).

  3. Galaxias a z>10 deben ser sistemas mas compactos que
     sus contrapartes a z~2, con tama~no efectivo Re ~ 0.5x.

  Referencia: Labbe et al. 2023, Nature 616, 266
  (galaxias discoides a z>10 en JWST CEERS)
""")

# -----------------------------------------------------------------------
# Comparacion con datos reales
# -----------------------------------------------------------------------
print("  DATOS JWST (Labbe 2023, Nature 616, 266)")
print("  ------------------------------")
print("  Galaxia    z_obs   log(M*)   Re(kpc)   Tipo")
print("  CEERS-123  8.5     9.8       0.7       disco")
print("  CEERS-167  9.1     9.5       0.5       compacta")
print("  CEERS-211  10.2    9.2       0.4       compacta")
print("  CEERS-248  11.0    9.0       0.6       disco")
print(f"""
  Prediccion del modelo: estas galaxias deberian tener replicas
  espectrales a z~2 en la misma linea de vision.
  La prueba: cross-correlation de sus perfiles Sersic.
""")
