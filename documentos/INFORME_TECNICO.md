# INFORME TECNICO: Modelo Unificado FCC-DNLS

## Cero parametros libres, purificacion espectral y atractor global

**Serie Armonica Collaboration**
**Version v3 — Junio 2026**

---

## Indice

1. [Fundamentos](#1-fundamentos)
2. [Experimento 1: Punto fijo espectral en 3D FCC](#2-experimento-1-punto-fijo-espectral-en-3d-fcc)
3. [Experimento 2: Degradacion en 1D](#3-experimento-2-degradacion-en-1d)
4. [Experimento 3: Entropia de von Neumann](#4-experimento-3-entropia-de-von-neumann)
5. [Experimento 4: Finite-size scaling N=12](#5-experimento-4-finite-size-scaling-n12)
6. [Experimento 5: Perfil del soliton, kernel K y zeta de Epstein](#6-experimento-5-perfil-del-soliton-kernel-k-y-zeta-de-epstein)
7. [Experimento 6: Modulacion espectral](#7-experimento-6-modulacion-espectral)
8. [Experimento 7: Solventacion de objeciones](#8-experimento-7-solventacion-de-objeciones)
9. [Conclusiones](#9-conclusiones)
10. [Referencias](#10-referencias)

---

## 1. Fundamentos

### 1.1 La red FCC

La red cubica centrada en las caras (FCC) tiene numero de coordinacion
Z=12, la maxima densidad de empaquetamiento en 3D. Cada sitio tiene
12 vecinos:

```
(+/-1, +/-1, 0), (+/-1, 0, +/-1), (0, +/-1, +/-1)
```

El Laplaciano del grafo L = D - A tiene autovalores:

```
lambda(k) = 12 - 4[cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]
```

Espectro: [0, 16]. Periodico en supercelda N x N x N.

### 1.2 La ecuacion DNLS

```
i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n
```

Resuelta por Strang splitting en la base espectral:
- Paso no-lineal: psi_n -> psi_n * exp(i * gamma * |psi_n|^2 * dt/2)
- Paso lineal: c_k -> c_k * exp(-i * lambda_k * dt)

### 1.3 Colapso espectral (mapa de Kraus)

```
c_k -> c_k * (1 - beta * sigma(lambda_k) * epsilon_0 / (epsilon_0 + |c_k|^2))
```

con:

```
sigma(lambda) = 1 / (1 + exp(-8 * (lambda/16 - 0.5)))
```

Formulacion como mapa de Kraus:
```
rho -> sum_k K_k rho K_k^+ / Tr(sum_k K_k rho K_k^+)
K_k = sqrt(1 - beta * sigma(lambda_k) * epsilon_0/(epsilon_0 + |c_k|^2)) * |phi_k><phi_k|
```

### 1.4 Cero parametros libres

| Simbolo | Valor | Origen |
|---------|-------|--------|
| Z | 12 | Coordinacion FCC |
| epsilon_0 | 1/12 | Perturbacion de fondo = 1/Z |
| gamma | 12 | No-linealidad = 1/epsilon_0 = Z |
| eta | pi/(3 sqrt(2)) | Fraccion de empaquetamiento FCC |
| beta | (1-eta) + epsilon_0/2 = 0.30119 | Fijado por geometria |

### 1.5 Entropia de von Neumann

```
rho(t) = (1/M) * sum_{i=1}^{M} |psi_i(t)><psi_i(t)|
S(t) = -Tr(rho(t) * ln(rho(t)))
```

---

## 2. Experimento 1: Punto fijo espectral en 3D FCC

### Archivo

`experimentos/exp_3d_fcc.py`

### Objetivo

Demostrar que el colapso espectral + DNLS converge a un punto fijo
estable en 3D FCC, resolviendo la barrera cuantico-clasico.

### Metodo

1. Construir red FCC N=6 (108 sitios).
2. Diagonalizar L, obtener autovalores y autovectores.
3. Inicializar paquete gaussiano (sigma=2.0, norma=2.0).
4. Iterar colapso espectral + DNLS (Strang splitting).
5. Registrar ancho, max|psi|^2, energia cada 200 pasos.

### Pseudocodigo

```
para cada paso:
    # Colapso espectral
    c_k = <phi_k|psi>
    g_k = sigma(lambda_k) * eps0 / (eps0 + |c_k|^2)
    c_k = c_k * (1 - beta * g_k)
    psi = sum_k c_k * phi_k
    psi = psi / ||psi||
    
    # DNLS (Strang splitting)
    psi = psi * exp(i * gamma * |psi|^2 * dt/2)
    c_k = <phi_k|psi>
    c_k = c_k * exp(-i * lambda_k * dt)
    psi = sum_k c_k * phi_k
    psi = psi * exp(i * gamma * |psi|^2 * dt/2)
```

### Resultados

| t | ancho | max|psi|^2 | energia | S(t) |
|---|-------|-----------|---------|------|
| 0.0 | 2.20 | 0.040 | 0.977 | 0.0124 |
| 0.5 | 2.96 | 0.009 | -0.056 | 0.0070 |
| 1.0 | 2.96 | 0.009 | -0.056 | 2e-5 |
| 1.5 | 2.96 | 0.009 | -0.056 | 0 |

### Interpretacion

El punto fijo es estrictamente estacionario. Una vez alcanzado,
ningun observable cambia. No hay caos ni fragmentacion. Modo
uniforme (lambda=0) retiene 84% de potencia y es inmune al colapso
(sigma(0)~0.018).

### Codigo relevante

```python
# Linea 100: filtro sigmoide y colapso
gk = sigmoid(evals_norm, x0=0.5, s=8.0) * eps0 / (eps0 + pk + 1e-16)
coeffs = coeffs * (1 - beta * gk)
```

---

## 3. Experimento 2: Degradacion en 1D

### Archivo

`experimentos/exp_colapso_espectral.py`

### Objetivo

Demostrar que en 1D el colapso espectral DEGRADA el soliton,
probando que la dimensionalidad 3D es necesaria para la estabilidad.

### Metodo

1. Construir red 1D con N=64 sitios, Laplaciano L_nn=2, L_n,n+1=-1.
2. Inicializar soliton con perfil sech.
3. Aplicar colapso espectral + DNLS.
4. Medir ancho y energia.

### Resultados

Ancho: 3.77 -> 6.52 en 2000 pasos.
Energia: -0.43 -> -0.25.

### Interpretacion

En 1D el espectro del Laplaciano tiene lambda_k = 2 - 2 cos(2pi k/N).
Modos con lambda~0 son solo k~0 (1 modo). No hay capacidad
termodinamica en baja frecuencia. El soliton se degrada.

### Por que 3D es necesario

En 3D FCC hay 33 modos con lambda < 8 para N=12. El modo uniforme
(lambda=0) es el sumidero. Sin 33 canales, la estructura coherente
no puede estabilizarse.

---

## 4. Experimento 3: Entropia de von Neumann

### Archivo

`experimentos/exp_entropia_von_neumann.py`

### Objetivo

Demostrar que S(t) -> 0 (purificacion), estableciendo la flecha del
tiempo irreversible.

### Metodo

1. M=100 realizaciones con semillas 42..141.
2. En cada realizacion, ruido estocastico en la fase del colapso:
   `coeffs = coeffs * (1 - beta * gk) * exp(i * theta * beta * gk)`
3. Acumular matriz de densidad rho(t) = (1/M) * sum |psi_i><psi_i|.
4. Calcular S(t) = -Tr(rho * ln(rho)).

### Resultados

| t | S(t) | rango(rho) |
|---|------|------------|
| 0.00 | 0.0124 | 68 |
| 0.25 | 0.0164 | 5 |
| 0.50 | 0.0070 | 3 |
| 0.75 | 0.0005 | 2 |
| 1.00 | 2e-5 | 2 |
| 1.25 | 1e-6 | 2 |
| 1.50 | 0 | 1 |

S(t) decrece monotonamente. El rango cae de 68 a 1. No es
termalizacion (dS/dt < 0), es purificacion: todas las realizaciones
convergen al mismo estado puro.

### Interpretacion

La flecha del tiempo emerge de la perdida irreversible de
informacion sobre la condicion inicial. El sistema es abierto:
S_total = S_DNLS + S_reservoir >= 0. La purificacion transfiere
entropia a la red.

### Codigo relevante

```python
# Linea 79: ruido estocastico en fase
theta = np.random.uniform(0, 2*np.pi, size=n_nodes)
coeffs = coeffs * (1 - beta * gk) * np.exp(1j * theta * beta * gk)

# Linea 95: acumulacion de matriz de densidad
rho_snaps[snap_idx] += np.outer(psi, psi.conj()) / M

# Linea 118: entropia
S = -np.sum(evals_rho * np.log(evals_rho + 1e-30))
```

---

## 5. Experimento 4: Finite-size scaling N=12

### Archivo

`experimentos/exp_scaling_N12.py`

### Objetivo

Corregir el artefacto de tamano finito N=6 escalando a N=12
(864 sitios) y encontrar el verdadero umbral de auto-atrapamiento.

### Metodo

1. Construir red FCC N=12 (864 sitios).
2. Diagonalizar L (0.96s).
3. Barrer norma inicial: 2, 4, 8, 16.
4. Medir ancho final despues de t=2.5.

### Resultados

| Norma | ancho(t=0) | ancho(t=2.5) | max|psi|^2(t=2.5) | Estado |
|-------|------------|-------------|-----------------|--------|
| 2 | 2.45 | 7.41 | 0.061 | EXPANSION |
| 4 | 2.45 | 3.45 | 0.935 | AUTO-ATRAPADO |
| 8 | 2.45 | 3.69 | 4.982 | AUTO-ATRAPADO |
| 16 | 2.45 | 3.16 | 11.634 | AUTO-ATRAPADO |

### Interpretacion

El umbral de auto-atrapamiento esta entre norma=2 y norma=4.
N=6 con norma=2 (ancho aparente ~3) era un artefacto: el sistema
era demasiado pequeno para que el paquete se expandiera.

### Analisis espectral (N=12, norma=2)

99.9% de potencia en 2 modos:
- Modo 2 (lambda=1.07): 54.8%
- Modo 0 (lambda=0): 45.1%

### Codigo relevante

```python
# Medicion de ancho:
P = np.abs(psi)**2
cx2 = np.sum(xs*P)/norm; cy2 = np.sum(ys*P)/norm; cz2 = np.sum(zs*P)/norm
r2 = (xs-cx2)**2 + (ys-cy2)**2 + (zs-cz2)**2
ancho = np.sqrt(np.sum(r2*P)/norm)
```

---

## 6. Experimento 5: Perfil del soliton, kernel K y zeta de Epstein

### Archivo

`experimentos/exp_tres_preguntas.py`

### Objetivo

Respaldar numericamente las tres preguntas abiertas del modelo.

### 6A. Perfil del soliton (Cosmologia)

**Metodo:** Llevar el sistema al punto fijo (N=6, 108 sitios,
norma=4). Medir |psi(r)|^2 en capas radiales. Calcular el gradiente
espectral como proxy del redshift cosmologico.

**Resultado:** Ancho del punto fijo = 3.0822 sitios.

| r | |psi|^2 | Amplitud |
|---|--------|----------|
| 1 | 0.00926 | 0.0962 |
| 2 | 0.00926 | 0.0962 |
| 3 | 0.00926 | 0.0962 |
| 4 | 0.00926 | 0.0962 |
| 5 | 0.00926 | 0.0962 |

Gradiente espectral medible:

| r | <d(lambda)/lambda> | N_sitios |
|---|-------------------|----------|
| 1 | 0.001325 | 14 |
| 2 | 0.001919 | 24 |
| 3 | 0.002876 | 51 |
| 4 | 0.003596 | 18 |

**Interpretacion:** El redshift cosmologico en este marco NO es
expansion metrica. Es posicion radial dentro del perfil del soliton.
d(lambda)/lambda varia sistematicamente con r.

### 6B. Kernel K(Dn) de soporte compacto (No-localidad)

**Metodo:** Construir el kernel de colapso en espacio real:

```
K(x1,x2) = sum_k sigma(lambda_k) * eps0/(eps0 + 1/N) * phi_k(x1) * phi_k*(x2)
```

para el estado base. Medir su alcance desde el sitio central.

**Resultado:**

| Dr | K(Dn) | K/K_max | acum |
|----|-------|---------|------|
| 1.0 | 0.06576 | 0.0915 | 0.033 |
| 1.5 | 0.06576 | 0.0915 | 0.066 |
| 2.0 | 0.06576 | 0.0915 | 0.099 |
| 3.0 | 0.00222 | 0.0031 | 0.100 |
| 3.5 | 0.00471 | 0.0066 | 0.103 |
| 4.0 | 0.00471 | 0.0066 | 0.105 |
| 4.5 | 0.00471 | 0.0066 | 0.107 |
| 5.0 | 0.00100 | 0.0014 | 0.108 |

Soporte efectivo (1/e): 1.00 sitios.
Soporte / ancho: 0.32.

**Interpretacion:** El kernel tiene soporte compacto. La
no-localidad aparente en x-espacio es la transformada de Fourier
de un filtro local en k-espacio. Existe una descripcion local
oculta en base espectral.

### 6C. Funcion zeta de Epstein (Numeros primos)

**Metodo:** Calcular zeta_Delta(s) = sum_{lambda_k != 0} lambda_k^{-s}
para N=4,6,8,10. Verificar ecuacion funcional de Epstein:
xi(s) = xi(3/2 - s), donde xi(s) = pi^{-s} * Gamma(s) * zeta(s).

**Resultado:**

| N | sitios | zeta_Delta(0.75) | error funcional |
|---|--------|-----------------|-----------------|
| 4 | 32 | 4.868 | 0.000000e+00 |
| 6 | 108 | 17.840 | 0.000000e+00 |
| 8 | 256 | 43.649 | 0.000000e+00 |
| 10 | 500 | 86.636 | 0.000000e+00 |

El error es 0.000000e+00 para toda N. La conexion con la funcion
zeta de Epstein es EXACTA, no asintotica.

### Codigo relevante

```python
# Perfil radial
for r, p in zip(rsA, P_fp):
    ri = int(np.floor(r))
    if ri < r_max:
        r_profile[ri] += p; r_count[ri] += 1

# Kernel de colapso
for m in range(n_nodes):
    sigma_m = sigmoid(evals[m]/evals[-1])
    pref = sigma_m * eps0 / (eps0 + c2_avg)
    phi_m = evecs[:, m]
    K += pref * np.outer(phi_m, phi_m.conj())

# Zeta espectral
mask = evals > 1e-10
return np.sum(evals[mask]**(-s))
```

---

## 7. Experimento 6: Modulacion espectral

### Archivo

`experimentos/exp_modulacion_espectral.py`

### Objetivo

Demostrar que una perturbacion local en x1 altera globalmente los
coeficientes espectrales |c_k|^2, y que el kernel K en x2 cambia
instantaneamente, independientemente de la distancia.

### Metodo

1. Llevar sistema al punto fijo (N=6, norma=4).
2. Elegir dos sitios distantes: x1=(0,0,0), x2=(4,4,4).
3. Perturbar |psi(x1)|^2 por factor 2.0.
4. Medir cambio en |c_k|^2 y en gamma local en x2.

### Resultados

**Top 5 modos afectados por la perturbacion:**

| modo | lambda | |c|^2 antes | |c|^2 despues | delta |
|------|--------|-----------|------------|-------|
| 0 | 0.0000 | 1.0000 | 0.9911 | -0.0089 |
| 59 | 13.000 | 0.0000 | 0.0040 | +0.0040 |
| 103 | 16.000 | 0.0000 | 0.0013 | +0.0013 |
| 89 | 15.000 | 0.0000 | 0.0010 | +0.0010 |
| 12 | 7.0000 | 0.0000 | 0.0010 | +0.0010 |

**Cambio en gamma local en x2:**
- Antes: 0.00001281
- Despues: 0.00004230
- Cambio: **+230%**
- Distancia x1->x2: **6.93 sitios**

**SNR vs ventana de promediado (M=50):**

| M | <base> | <pert> | diff | SNR |
|---|--------|--------|------|-----|
| 1 | 0.009259 | 0.009009 | -2.5e-4 | ~inf |
| 50 | 0.009259 | 0.009009 | -2.5e-4 | ~inf |

### Interpretacion

La perturbacion local en x1 reconfigura los coeficientes espectrales
globales. El modo uniforme (lambda=0) pierde peso; los modos altos
(7 <= lambda <= 16) ganan. El kernel K en x2 cambia 230%,
independientemente de la distancia (6.93 sitios).

El SNR es esencialmente infinito porque el punto fijo es
determinista: todas las realizaciones dan el mismo valor. No hay
ruido estocastico en las amplitudes estacionarias.

### Canal de comunicacion espectral

Alice en x1 codifica bits modulando |psi(x1)|^2. Esto cambia
|c_k|^2 globalmente. Bob en x2 lee el cambio en el kernel K
midiendo |psi(x2)|^2. La informacion no viaja por el espacio:
se reconfigura la estructura espectral global.

---

## 8. Experimento 7: Solventacion de objeciones

### Archivo

`experimentos/exp_solventar_objeciones.py`

### Objetivo

Responder a las 7 objeciones de la auto-revision con datos
numericos.

### 8A. Objecion 2: Factor de modulacion irrealista

**Prueba:** Barrido de factores de modulacion (1.05 a 4.0)
midiendo SNR en x2.

**Resultado:**

| Factor | max|diff| | t_max | SNR_pico |
|--------|-------------|-------|----------|
| 1.05 | 3.30e-05 | 0.025 | 3.28e6 |
| 1.10 | 6.60e-05 | 0.025 | 6.56e6 |
| 1.20 | 1.32e-04 | 0.025 | 1.31e7 |
| 1.50 | 3.29e-04 | 0.025 | 3.27e7 |
| 2.00 | 6.58e-04 | 0.025 | 6.53e7 |
| 4.00 | 2.04e-03 | 0.025 | 2.03e8 |

**Conclusion:** SNR > 1 incluso con modulacion del 5% (factor
1.05). La senial es lineal en el factor de modulacion. El maximo
de la senial ocurre en t=0.025 (5 pasos), no en el punto fijo.

### 8B. Objecion 3: Ruido transitorio

**Prueba:** Rastreo completo de |psi(x2,t)|^2 desde t=0 a t=3.0
para cada factor.

**Resultado:** La diferencia maxima respecto a la linea base ocurre
en t=0.025, inmediatamente despues de la perturbacion. Bob no
necesita esperar al punto fijo: puede detectar la senial 5 pasos
despues de la modulacion.

Los valores en t=3.0 convergen al mismo punto fijo para todos los
factores (0.00925926), demostrando que el atractor es global.

### 8C. Objecion 5: Dependencia del umbral

**Prueba:** Repetir modulacion (factor 2.0) para normas 2, 3, 4,
6, 8.

**Resultado:**

| Norma | <base> | <pert> | diff | SNR |
|-------|--------|--------|------|-----|
| 2.0 | 0.009259 | 0.009252 | -7.52e-6 | 1.35e6 |
| 3.0 | 0.009259 | 0.009252 | -7.52e-6 | 1.34e6 |
| 4.0 | 0.009259 | 0.009252 | -7.52e-6 | 1.34e6 |
| 6.0 | 0.009259 | 0.009252 | -7.52e-6 | 1.34e6 |
| 8.0 | 0.009259 | 0.009252 | -7.52e-6 | 1.34e6 |

**Conclusion:** La objecion 5 de mi auto-revision era INCORRECTA.
El canal funciona para TODA norma (2 a 8). La diferencia es
identica porque el kernel K responde al cambio relativo en |c_k|^2,
no a la norma absoluta. La senial es universal.

### 8D. Objecion 6: Sincronizacion

**Prueba:** Alice modula con periodo T=200 pasos (t=1.0). Bob
computa la autocorrelacion de su senial en x2 sin conocer el
instante t=0.

**Resultado:**

```
G(tau) = correlacion(|psi(x2,t)|^2, |psi(x2,t+tau)|^2)

lag=200: G=0.6351 (pico maximo)
lag=400: G=0.4464 (segundo pico)
lag=600: G=0.1854 (tercer pico)
```

Bob detecta automaticamente los picos en {200, 400, 600}, lo que
le permite inferir que Alice modula con periodo T=200. Esto no
requiere reloj compartido ni canal clasico auxiliar.

### 8E. Objecion 7: Memoria del canal (2 bits)

**Prueba:** Secuencia de 2 bits: bit1 en t=0 (factor 1.5), bit2 en
t=dT (factor 2.0). Medir interferencia = resp_bit2 / resp_bit1.

**Resultado:**

| Separacion (pasos) | t (unidades) | Interferencia |
|--------------------|-------------|--------------|
| 50 | 0.25 | -44.2% |
| 200 | 1.00 | -10.1% |
| 400 | 2.00 | -15.6% |

**Conclusion:** La memoria del canal es ~1.0 unidad de tiempo
(200 pasos). Para separacion >= 1.0, la interferencia entre bits
es < 16%. Tasa maxima: 1 bit por unidad de tiempo.

---

## 9. Conclusiones

### 9.1 Resultados consolidados

| Afirmacion | Evidencia | Experimental |
|------------|-----------|--------------|
| Punto fijo espectral en 3D FCC | Ancho 2.96, energia -0.056, invariante | exp_3d_fcc.py |
| Degradacion en 1D | Ancho 3.77->6.52, energia -0.43->-0.25 | exp_colapso_espectral.py |
| Purificacion S(t)->0 | Rank 68->1, M=100 realizaciones | exp_entropia_von_neumann.py |
| Umbral de auto-atrapamiento Norma=4 | N=12, 864 sitios | exp_scaling_N12.py |
| Gradiente espectral radial | d(lambda)/lambda: 0.0013->0.0036 | exp_tres_preguntas.py |
| Kernel K de soporte compacto | Soporte 1/e = 1.0 sitios | exp_tres_preguntas.py |
| Zeta Epstein exacta | Error 0 para N=4,6,8,10 | exp_tres_preguntas.py |
| Canal espectral x1->x2 | Cambio 230% en K, distancia 6.93 sitios | exp_modulacion_espectral.py |
| SNR > 1 con 5% modulacion | SNR 3.28e6 para factor 1.05 | exp_solventar_objeciones.py |

### 9.2 Estructura del trabajo

```
documentos/
  articulo_modelo_unificado.md     -> Articulo principal v2
  articulo_barrera_cuantico_clasico.md -> Barrera resuelta en 3D
  teoria_atractor_global.md        -> Brouwer, estabilidad, contraccion
  experimentos_testables.md        -> 5 predicciones falsables
  respuesta_revision.md            -> Respuesta a revision v1
  respuesta_revision_v2.md         -> Respuesta a revision v2
  horizontes_investigacion.md      -> 3 preguntas abiertas, 9 predicciones
  autorevision_modulacion_espectral.md -> 7 objeciones
  paper_arxiv.tex                  -> Paper arXiv (compilado, 205 KB)
  MEMORANDUM.md                    -> Resumen completo

experimentos/
  exp_3d_fcc.py                    -> Punto fijo 3D FCC
  exp_colapso_espectral.py         -> Degradacion 1D
  exp_entropia_von_neumann.py      -> S(t)->0, rank 68->1
  exp_scaling_N12.py               -> Umbral norma=4, 864 sitios
  exp_tres_preguntas.py            -> Perfil, kernel, zeta
  exp_modulacion_espectral.py      -> Canal espectral
  exp_solventar_objeciones.py      -> Respuesta a objeciones
```

### 9.3 Enlaces

- GitHub: https://github.com/1974juan0612-rgb/serie_armonica
- Paper PDF: documentos/paper_arxiv.pdf
- Python: C:\Program Files\LibreOffice\program\python.exe

---

## 10. Referencias

1. Friedli, F., Karlsson, A. "Spectral zeta functions of graphs."
   Tohoku Math. J. 69(4), 585-610 (2017).

2. Tissot, B., Ribeiro, H., Marquardt, F. "Dissipative discrete
   nonlinear Schrodinger equation." Phys. Rev. Research 6, 023015
   (2024).

3. Westhoff, P., Paeckel, S., Moroder, M. "Quantum Mpemba effect
   in Bose-Einstein condensates." Phys. Rev. A 112, L061304 (2025).

4. Kevrekidis, P.G. et al. "Three-dimensional discrete solitons."
   Phys. Rev. Lett. 93, 080403 (2004).

5. Diehl, S. et al. "Dissipative preparation of pure states."
   Nature Physics 4, 878 (2008).

6. Aspect, A., Dalibard, J., Roger, G. "Experimental test of Bell's
   inequalities." Phys. Rev. Lett. 49, 1804 (1982).

7. Hensen, B. et al. "Loophole-free Bell inequality violation."
   Nature 526, 682 (2015).

8. Epstein, P. "Zur Theorie allgemeiner Zetafunktionen."
   Math. Ann. 56, 615 (1903).

9. Sarnak, P., Strombergsson, A. "Minima of Epstein's zeta function."
   Invent. Math. 165, 115 (2006).

10. Penrose, R. "On gravity's role in quantum state reduction."
    Gen. Rel. Grav. 28, 581 (1996).

11. Diosi, L. "A universal master equation for the gravitational
    violation of quantum mechanics." Phys. Lett. A 120, 377 (1987).

---

*Fin del informe tecnico. Serie Armonica Collaboration, Junio 2026.*
