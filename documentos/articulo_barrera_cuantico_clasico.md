# Emergencia de la barrera cuantico-clasico en un modelo de red discreta con auto-interaccion

**Autor:** Serie Armonica Collaboration

**Fecha:** Junio 2026 (v2: incorporacion de separacion espectral 3D)

## Resumen

Presentamos un modelo de espacio-tiempo discreto basado en una red FCC con ecuacion
de Schrodinger no-lineal discreta (DNLS) y un mecanismo de colapso auto-regulado en el
espacio espectral. El modelo reproduce la barrera entre el comportamiento cuantico
(colapsos probabilisticos en alta frecuencia) y el clasico (solitones estables en baja
frecuencia) sin postular la dicotomia, sino viendola emerger de la dinamica. En 1D,
los colapsos y la no-linealidad compiten destructivamente generando caos o sobre-enfoque.
En 3D FCC, la separacion espectral revela un punto fijo estable donde ambos mecanismos
coexisten: el modo uniforme (lambda=0) actua como sumidero termodinamico protegiendo la
estructura clasica mientras los modos de alta frecuencia colapsan como ruido cuantico.
La flecha del tiempo emerge como purificacion irreversible del ensemble (entropia de
von Neumann S(t) -> 0). Este resultado sugiere que el problema de la medida cuantica
no es un artefacto filosofico sino una propiedad emergente de sistemas con auto-interaccion
en redes discretas, y que su resolucion requiere dimensionalidad >= 3.

## 1. Introduccion

El problema de la medida en mecanica cuantica -- la transicion entre la superposicion
de estados y un resultado clasico definido -- ha sido discutido durante casi un siglo
sin resolucion definitiva. Aqui mostramos que este problema emerge naturalmente en un
modelo minimalista: una red FCC con ecuacion de Schrodinger no-lineal discreta y un
mecanismo de colapso espectral auto-regulado. La dependencia dimensional del resultado
(1D: caos, 3D: punto fijo estable) sugiere que la solucion al problema de la medida
requiere un espacio con al menos tres dimensiones.

## 2. El modelo

### 2.1 Evolucion unitaria (DNLS)

    i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n                     (1)

con L = D - A el Laplaciano de la red FCC (12 vecinos/sitio, autovalores en [0, 16]),
gamma = 12 = Z_FCC (determinado por la geometria, ver articulo companion).

### 2.2 Colapso espectral auto-regulado

En la base de autovectores {phi_k, lambda_k} del Laplaciano, expandimos
psi = sum_k c_k phi_k. El colapso actua selectivamente:

    c_k --> c_k * (1 - beta * sigma(lambda_k) * epsilon_0/(epsilon_0 + |c_k|^2))   (2)

con sigma(lambda) = 1/[1 + exp(-8(lambda/16 - 0.5))], epsilon_0 = 1/12, beta = 0.3.

Los modos con lambda << 8 (baja frecuencia, estructura) tienen sigma ~ 0 y no colapsan.
Los modos con lambda >> 8 (alta frecuencia, ruido) tienen sigma ~ 1 y colapsan
con tasa regulada por epsilon_0/(epsilon_0 + |c_k|^2).

### 2.3 Entropia de von Neumann

La matriz de densidad promediada sobre M realizaciones estocasticas:

    rho(t) = (1/M) sum_i |psi_i(t)><psi_i(t)|                                     (3)
    S(t) = -Tr(rho(t) ln rho(t))                                                   (4)

S(t) mide la pureza del ensemble. S = 0 indica un estado puro (todas las
realizaciones identicas). S > 0 indica mezcla estadistica.

## 3. Resultados

### 3.1 Soliton puro (sin colapsos)

La DNLS con gamma = 12 forma solitones estables desde pulsos localizados:

    Ancho inicial: 3.77
    Ancho en t=6:  2.74  (compresion del 27%)
    Energia: conservada a 10^-6

### 3.2 Colapso en espacio real (V1: modelo original)

Cuando el colapso actua en espacio real (sitio por sitio), dos regimenes:

**Colapso aleatorio** (psi_n -> ruido con prob. gamma_n = eps0/(eps0+|psi_n|^2)):

    Ancho oscila entre 7 y 22 sitios
    Energia fluctua: -0.78 < E < 0.65
    Sin convergencia a estado estable

**Colapso filtrante** (psi_n * (1 - beta*gamma_n)):

    Ancho = 0.0075, max|psi|^2 = 0.9999
    Sobre-enfoque: toda la probabilidad en 1 sitio

Ambos regimenes son insatisfactorios: caos o sobre-enfoque. La barrera parece infranqueable.

### 3.3 Colapso espectral en 3D FCC: punto fijo estable

Cuando el colapso actua en la base espectral (ecuacion 2) en 3D FCC (N=6, 108 sitios):

| t | ancho | max|psi|^2 | energia |
|---|-------|-----------|---------|
| 0.0 | 2.20 | 0.040 | 0.977 |
| 1.0 | 2.96 | 0.009 | -0.056 |
| 5.0 | 2.96 | 0.009 | -0.056 |

PUNTO FIJO ESTABLE en t ~ 1.5. Colapso y DNLS coexisten sin caos ni sobre-enfoque.

### 3.4 Dependencia dimensional: 1D vs 3D

**1D (N=64)**: el soliton se DEGRADA bajo colapso espectral.
Ancho crece de 3.77 a 6.52 en 2000 pasos. No hay punto fijo.
Razon: el soliton 1D (perfil sech) distribuye su potencia en muchos modos de Fourier,
y ninguno alcanza la amplitud necesaria para protegerse del colapso.

**3D FCC (N=6, 108 sitios)**: punto fijo ESTABLE.
Razon: el modo uniforme (lambda=0, sigma(0)~0.018) retiene el 84% de la potencia
y es practicamente inmune al colapso. El modo uniforme actua como sumidero termodinamico.

### 3.5 Finite-size scaling: N=12 (864 sitios)

Barrido de la norma del paquete gaussiano inicial (sigma=2.0):

| Norma | ancho(t=2.5) | Estado |
|-------|-------------|--------|
| 2 | 7.41 | EXPANSION |
| 4 | 3.45 | AUTO-ATRAPADO |

Umbral de auto-atrapamiento en norma=4. Por debajo, la dispersion domina.
Por encima, la no-linealidad forma estructura coherente.

### 3.6 Estructura espectral del punto fijo

N=12, norma=2: el punto fijo tiene el 99.9% de la potencia en los dos modos
mas bajos (lambda=0 y lambda=1.07). Solo 2 modos significativos de 864 posibles.
El colapso espectral "elige" los modos de minima frecuencia como nucleo clasico.

### 3.7 Entropia de von Neumann

M=100 realizaciones, 3D FCC (N=6):

| t | S(t) | rango(rho) |
|---|------|-----------|
| 0.00 | 0.0124 | 68 |
| 0.50 | 0.0070 | 3 |
| 1.00 | 2e-5 | 2 |
| 1.50 | 0 | 1 |

S(t) -> 0 monotonicamente. El ensemble se PURIFICA: todas las realizaciones
convergen al mismo estado puro. Este es el analogo dinamico del colapso de
la funcion de onda en la medicion cuantica.

## 4. Interpretacion: la barrera emerge, la separacion espectral la resuelve

### 4.1 Exclusion mutua en espacio real

Colapso y no-linealidad en espacio real son mutuamente excluyentes:

    Colapso aleatorio + DNLS = caos
    Filtro + DNLS = sobre-enfoque
    Solo DNLS = soliton
    Solo colapso = ruido

### 4.2 Coexistencia en espacio espectral 3D

La separacion espectral permite la coexistencia:

    Colapso en alta frecuencia + DNLS en baja frecuencia = PUNTO FIJO

La clave es que el Laplaciano FCC tiene un modo uniforme (lambda=0) que es
inmune al colapso espectral (sigma(0) ~ 0). Este modo ancla la estructura
clasica mientras los modos de alta frecuencia fluctuan como ruido cuantico.

### 4.3 La flecha del tiempo es purificacion

S(t) -> 0 no es termalizacion (dS/dt > 0) sino purificacion (dS/dt < 0).
La direccion del tiempo esta definida por la perdida irreversible de mezcla
estadistica: desde un ensemble heterogeneo (68 autoestados ocupados) hacia
un estado puro unico (rango=1). Esto es exactamente lo que ocurre en una
medicion cuantica.

## 5. Proxima direccion

### 5.1 Formalizacion del atractor global

El punto fijo espectral es un atractor en el espacio de Hilbert (no solo
en el espacio de configuraciones). Demostrar su unicidad y estabilidad
Lyapunov es el siguiente paso teorico.

### 5.2 Continuo no-lineal

La conexion entre el colapso espectral discreto y la ecuacion de Schrodinger
no-lineal continua (cubic NLS con ruido disipativo) es una direccion natural
para tender puentes con la literatura establecida.

### 5.3 Implicaciones para la medicion cuantica

Si la separacion de escalas via el Laplaciano del espacio fisico resuelve
la barrera cuantico-clasico, entonces el problema de la medida no requiere
nueva fisica (colapso objetivo, gravedad cuantica, etc.) sino la comprension
de que la no-linealidad y el colapso operan en subespacios espectrales
complementarios del Laplaciano del espacio.

## 6. Conclusion

Un modelo minimalista de red FCC con DNLS y colapso espectral auto-regulado
reproduce la barrera cuantico-clasico y la resuelve en 3D mediante separacion
de escalas en el espectro del Laplaciano. La flecha del tiempo emerge como
purificacion irreversible del ensemble. Esto sugiere que el problema de la
medida cuantica no es externo a la teoria, sino una consecuencia inevitable
de la geometria discreta del espacio, y que su resolucion requiere la
dimensionalidad 3D del mundo real.

## Referencias

- Penrose, R. "On Gravity's Role in Quantum State Reduction." Gen. Rel. Grav. 28,
  581-600 (1996).
- Diosi, L. "A Universal Master Equation for the Gravitational Violation of Quantum
  Mechanics." Phys. Lett. A 120, 377-381 (1987).
- Friedli, F., Karlsson, A. "Spectral zeta functions of graphs." Tohoku Math. J.
  69(4), 585-610 (2017).
- Tsironis, G. "The Discrete Nonlinear Schrodinger Equation." Springer (2025).
- Connes, A., Consani, C. "Zeta Spectral Triples." arXiv:2511.22755 (2025).
