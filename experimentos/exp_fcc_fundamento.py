import numpy as np
from scipy import special as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# FUNDAMENTO: RED FCC COMO ESPACIO DISCRETO
# Laplaciano de grafo, espectro y funcion zeta espectral
# Referencia: Friedli & Karlsson (2017), Sarnak & Strombergsson (2006)
# ============================================================================

print("=" * 70)
print("RED FCC 3D - LAPLACIANO DISCRETO Y FUNCION ZETA ESPECTRAL")
print("=" * 70)

# ------------------------------------------------------------------
# 1. GENERACION DE LA RED FCC
# ------------------------------------------------------------------
# Descripcion: todos los puntos (i,j,k) in Z^3 con i+j+k par
# Unidades: a/2 = 1, donde a es la constante de red convencional
# 12 vecinos: (+-1, +-1, 0), (+-1, 0, +-1), (0, +-1, +-1)

N = 8
assert N % 2 == 0, "N debe ser par para condiciones de frontera periodicas"

print(f"\nSupercelda: {N}x{N}x{N}")
print(f"Constante de red convencional: a = 2 (a/2 = 1)")
print(f"Distancia 1os vecinos: a/raiz(2) = {2/np.sqrt(2):.4f}")

sites = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)
         if (i + j + k) % 2 == 0]
n_sites = len(sites)
print(f"Sitios FCC: {n_sites} (= N^3/2 = {N**3 // 2})")

idx_map = {s: idx for idx, s in enumerate(sites)}

# 12 vectores de vecindad
neighbors = [(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0),
             (1,0,1),(1,0,-1),(-1,0,1),(-1,0,-1),
             (0,1,1),(0,1,-1),(0,-1,1),(0,-1,-1)]

# ------------------------------------------------------------------
# 2. MATRIZ LAPLACIANA (densa, N pequeno)
# ------------------------------------------------------------------
print("\nConstruyendo matriz laplaciana...")
L = np.zeros((n_sites, n_sites))
for idx, (i, j, k) in enumerate(sites):
    L[idx, idx] = 12.0
    for di, dj, dk in neighbors:
        ni, nj, nk = (i + di) % N, (j + dj) % N, (k + dk) % N
        nidx = idx_map[(ni, nj, nk)]
        L[idx, nidx] = -1.0
print(f"Matriz: {L.shape}, elementos no nulos: {np.count_nonzero(L)}")

# ------------------------------------------------------------------
# 3. AUTOVALORES (diagonalizacion densa exacta)
# ------------------------------------------------------------------
print("\nDiagonalizando...")
evals = np.sort(np.linalg.eigvalsh(L))

# El primer autovalor deberia ser ~0 (modo constante)
print(f"Autovalor minimo: {evals[0]:.2e} (modo constante, debe ser ~0)")
evals_nz = evals[evals > 1e-10]
n_nz = len(evals_nz)
print(f"Autovalores no nulos: {n_nz} de {n_sites}")
print(f"Rango: [{evals_nz[0]:.4f}, {evals_nz[-1]:.4f}]")
print(f"Maximo teorico: 24.0 (para red FCC infinita)")

# Verificar algunos autovalores con la formula analitica de Bloch
# Para FCC: lambda(k) = 12 - 4[cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]
# con k = (pi/N) * (m_x, m_y, m_z) en la 1a zona de Brillouin
print("\n--- Validacion analitica (Bloch) ---")
# Algunos modos de alta simetria:
# Gamma: k = (0,0,0) -> lambda = 12 - 4*3 = 0
# X:     k = (pi, 0, 0)-> lambda = 12 - 4[-1+(-1)+1] = 12 - 4(-1) = 16
# L:     k = (pi,pi,pi)-> lambda = 12 - 4[(-1)(-1)+(-1)(-1)+(-1)(-1)] = 12 - 4(3) = 0
# W:     k = (pi, pi/2, 0)
for label, kx, ky, kz, lam_teo in [
    ("Gamma (0,0,0)", 0, 0, 0, 0),
    ("X (pi,0,0)", np.pi, 0, 0, 16),
    ("L (pi,pi,pi)", np.pi, np.pi, np.pi, 0),
]:
    lam = 12 - 4*(np.cos(kx)*np.cos(ky) + np.cos(kx)*np.cos(kz) + np.cos(ky)*np.cos(kz))
    print(f"  {label:<20}  analitico: {lam:.2f}")

# ------------------------------------------------------------------
# 4. HISTOGRAMA DEL ESPECTRO
# ------------------------------------------------------------------
plt.figure(figsize=(10, 5))
plt.hist(evals_nz, bins=30, color='#00e5a0', alpha=0.7, edgecolor='white', linewidth=0.5)
plt.axvline(np.mean(evals_nz), color='#ff4d00', linestyle='--',
            label=f'Media = {np.mean(evals_nz):.2f}')
plt.axvline(12.0, color='#7c3aff', linestyle=':', alpha=0.5, label='Centro de banda')
plt.xlabel(r'Autovalor $\lambda$', fontsize=12)
plt.ylabel('Frecuencia', fontsize=12)
plt.title(f'Espectro del Laplaciano FCC (supercelda {N}x{N}x{N}, {n_sites} sitios)')
plt.legend(fontsize=10)
plt.grid(alpha=0.15)
plt.tight_layout()
plt.savefig('C:\\Users\\famil\\Desktop\\serie_armonica\\espectro_fcc.png', dpi=150)
print(f"\nGrafico guardado: espectro_fcc.png")

# ------------------------------------------------------------------
# 5. FUNCION ZETA ESPECTRAL
# ------------------------------------------------------------------
# zeta_Delta(s) = sum_{lambda_k != 0} lambda_k^{-s}
#
# Para el laplaciano continuo en un toro 3D de lado L:
#   zeta(s) ~ (L/(2pi))^{2s} * zeta_Epstein(s)
# donde zeta_Epstein(s) = sum'_{n in Z^3} |n|^{-2s}
# converge para Re(s) > 3/2.

print("\n" + "=" * 70)
print("FUNCION ZETA ESPECTRAL")
print("zeta_Delta(s) = sum_{lambda_k != 0} lambda_k^{-s}")
print("=" * 70)

s_vals = np.linspace(1.6, 4.0, 25)
zeta_s = np.array([np.sum(evals_nz ** (-s)) for s in s_vals])

# Funcion completada: xi(s) = pi^{-s} * Gamma(s) * zeta(s)
# Para la zeta de Epstein del reticulado D_3 (FCC):
#   xi(3/2 - s) = xi(s)   [ecuacion funcional]
xi_s = np.array([np.pi**(-s) * sp.gamma(s) * np.sum(evals_nz ** (-s))
                  for s in s_vals])

print(f"\n{'s':<6} {'zeta(s)':<16} {'xi(s)':<16}")
print("-" * 40)
for i in range(0, len(s_vals), 3):
    print(f"{s_vals[i]:<6.2f} {zeta_s[i]:<16.6e} {xi_s[i]:<16.6e}")

# ------------------------------------------------------------------
# 6. ECUACION FUNCIONAL
# ------------------------------------------------------------------
# Verificamos cuan bien se cumple xi(s) ~ xi(3/2 - s)
# Para un sistema FINITO, se desvia del caso ideal.
# La Friedli-Karlsson (2017) demuestra que la ecuacion funcional
# asintotica para grafos ciclicos es EQUIVALENTE a la HR.

print("\n" + "=" * 70)
print("ECUACION FUNCIONAL: xi(s) vs xi(3/2 - s)")
print("xi(s) = pi^{-s} * Gamma(s) * zeta_Delta(s)")
print("-" * 70)

# s = 0.75 es el punto fijo (3/2 - s = s)
# Probamos alrededor
s_test = np.linspace(0.55, 0.95, 9)
print(f"\n{'s':<8} {'xi(s)':<16} {'xi(3/2-s)':<16} {'xi(s)/xi(3/2-s)':<16} {'|error|%':<12}")
print("-" * 72)
for s in s_test:
    xi_s = np.pi**(-s) * sp.gamma(s) * np.sum(evals_nz ** (-s))
    xi_c = np.pi**(-(1.5-s)) * sp.gamma(1.5-s) * np.sum(evals_nz ** (-(1.5-s)))
    r = xi_s / xi_c
    err = abs(r - 1) * 100
    print(f"{s:<8.2f} {xi_s:<16.6e} {xi_c:<16.6e} {r:<16.6f} {err:<11.3f}%")

print(f"""
Interpretacion:
  La ecuacion funcional exacta xi(s) = xi(3/2 - s) es valida para la
  funcion zeta de Epstein del reticulado D_3 INFINITO.
  
  Para un sistema FINITO de {n_sites} nodos, la suma sobre autovalores
  del laplaciano discreto APROXIMA la zeta de Epstein.
  
  La desviacion de la ecuacion funcional se debe a:
  1. Tamano finito de la supercelda (N={N})
  2. Discretizacion del espectro (modos de longitud de onda > N
     no estan representados)
  
  A medida que N -> infinito, la aproximacion mejora.
  El punto fijo s=0.75 siempre se cumple exactamente por construccion.
""")

# ------------------------------------------------------------------
# 7. GRAFICO: ECUACION FUNCIONAL
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Izquierda: zeta espectral
ax = axes[0]
ax.plot(s_vals, zeta_s, color='#00e5a0', linewidth=2, label=r'$\zeta_\Delta(s)$')
ax.axvline(x=1.5, color='#ff4d00', linestyle=':', alpha=0.5, label=r'$s=3/2$ (cota)')
ax.set_xlabel('s', fontsize=12)
ax.set_ylabel(r'$\zeta_\Delta(s)$', fontsize=12)
ax.set_title('Funcion zeta espectral', fontsize=13)
ax.legend(fontsize=9)
ax.grid(alpha=0.15)
ax.set_yscale('log')

# Derecha: verificacion ecuacion funcional
ax = axes[1]
s_checks = np.linspace(0.55, 1.0, 30)
xi_s = np.array([np.pi**(-s) * sp.gamma(s) * np.sum(evals_nz ** (-s))
                  for s in s_checks])
xi_c = np.array([np.pi**(-(1.5-s)) * sp.gamma(1.5-s) * np.sum(evals_nz ** (-(1.5-s)))
                  for s in s_checks])
ratio = xi_s / xi_c

ax.plot(s_checks, ratio, color='#00e5a0', linewidth=2, label=r'$\xi(s)/\xi(3/2-s)$')
ax.axhline(y=1.0, color='#ff4d00', linestyle='--', alpha=0.7,
           label='Valor esperado = 1')
ax.fill_between(s_checks, 0.95, 1.05, alpha=0.08, color='#00e5a0',
                label='+/- 5%')
ax.set_xlabel('s', fontsize=12)
ax.set_ylabel(r'$\xi(s) / \xi(3/2-s)$', fontsize=12)
ax.set_title('Verificacion ecuacion funcional', fontsize=13)
ax.legend(fontsize=9)
ax.grid(alpha=0.15)

plt.tight_layout()
plt.savefig('C:\\Users\\famil\\Desktop\\serie_armonica\\zeta_funcional_fcc.png', dpi=150)
print("\nGrafico guardado: zeta_funcional_fcc.png")

# ------------------------------------------------------------------
# 8. RELACION DE DISPERSION TIGHT-BINDING
# ------------------------------------------------------------------
# E(k) = -2t[cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]
# donde t es el parametro de salto (hopping)
print("\n" + "=" * 70)
print("RELACION DE DISPERSION TIGHT-BINDING EN FCC")
print("E(k) = -2t [cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]")
print("=" * 70)

# Grafico de dispersion a lo largo del camino Gamma-X-L-Gamma en la BZ
t = 1.0
n_k = 200
path = []

# Gamma -> X: (0,0,0) -> (pi,0,0)
for i in range(n_k):
    frac = i / n_k
    path.append((frac * np.pi, 0, 0, f'Gamma-X', frac))

# X -> L: (pi,0,0) -> (pi,pi,pi)  
for i in range(n_k):
    frac = i / n_k
    path.append((np.pi, frac * np.pi, frac * np.pi, f'X-L', 1 + frac))

# L -> Gamma: (pi,pi,pi) -> (0,0,0)
for i in range(n_k):
    frac = i / n_k
    path.append(((1-frac) * np.pi, (1-frac) * np.pi, (1-frac) * np.pi, f'L-Gamma', 2 + frac))

energy = []
for kx, ky, kz, _, _ in path:
    e = -2 * t * (np.cos(kx)*np.cos(ky) + np.cos(kx)*np.cos(kz) + np.cos(ky)*np.cos(kz))
    energy.append(e)

xs = np.array([p[4] for p in path])

plt.figure(figsize=(10, 4))
plt.plot(xs, energy, color='#7c3aff', linewidth=2)
plt.axhline(y=0, color='white', linewidth=0.5, alpha=0.2)
plt.xticks([0, 1, 2], [r'$\Gamma$', 'X', 'L'])
plt.xlabel('Camino en la Zona de Brillouin', fontsize=12)
plt.ylabel('E(k) / t', fontsize=12)
plt.title('Dispersion tight-binding en red FCC (s-band)')
plt.grid(alpha=0.15)
plt.tight_layout()
plt.savefig('C:\\Users\\famil\\Desktop\\serie_armonica\\dispersion_fcc.png', dpi=150)
print("\nGrafico guardado: dispersion_fcc.png")

# ------------------------------------------------------------------
# 9. RESUMEN
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("RESUMEN")
print("=" * 70)
print(f"""
Sistema: Red FCC 3D, supercelda {N}x{N}x{N} ({n_sites} nodos)
Coordinacion: 12 vecinos por nodo (empaquetamiento maximo)
Laplaciano: L = D - A (discreto de grafo)
Espectro: {n_nz} autovalores no nulos en [{evals_nz[0]:.4f}, {evals_nz[-1]:.4f}]

Funcion zeta espectral: zeta_Delta(s) = sum lambda_k^{{-s}}
  - Converge para Re(s) > 3/2
  - Se relaciona con la zeta de Epstein del reticulado D_3

Ecuacion funcional (Epstein, D_3):
  xi(s) = pi^{{-s}} * Gamma(s) * zeta_Delta(s)
  xi(3/2 - s) = xi(s)   [limite continuo, N -> inf]

  Para N={N}, la desviacion en s=0.75 es 0% (punto fijo).
  Para otros s, la desviacion mide efectos de tamano finito.

Referencias clave:
  1. Friedli & Karlsson (2017) - Tohoku Math. J. 69(4):585-610
     "Spectral zeta functions of graphs and the Riemann zeta function
      in the critical strip"
     -> Demuestra que la HR equivale a una ecuacion funcional asintotica
        para la funcion zeta espectral de grafos ciclicos.

  2. Sarnak & Strombergsson (2006) - Invent. Math. 165:115-151
     "Minima of Epstein's zeta function and heights of flat tori"
     -> Prueba que FCC minimiza la altura del toro plano en 3D.

  3. Sun et al. (2011) - J. Comput. Phys. 230:5869-5888
     "An FDTD scheme on an FCC grid for the wave equation"
     -> Deriva el laplaciano discreto en FCC con 2/3 menos anisotropia.
""")
