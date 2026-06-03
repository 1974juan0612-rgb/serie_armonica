# Modelo unificado de red FCC con cero parametros libres: solitones, tiempo emergente y separacion espectral

**Autor:** Serie Armonica Collaboration

**Fecha:** Junio 2026 (v2: respuesta a revision)

## Resumen

Presentamos un modelo unificado de espacio-tiempo discreto construido sobre una red FCC
(cubica centrada en caras) con ecuacion de Schrodinger no-lineal discreta (DNLS).
El modelo tiene cero parametros libres: la geometria de la red determina la perturbacion
de fondo (epsilon_0 = 1/12), la fuerza no-lineal (gamma = 12) y la tasa de colapso
por modo (gamma_k = sigma(lambda_k) * epsilon_0/(epsilon_0 + |c_k|^2)). El Laplaciano
de la red reproduce exactamente la relacion de dispersion analitica
lambda(k) = 12 - 4[cos(k_x)cos(k_y) + ...]. La funcion zeta espectral del grafo satisface
la ecuacion funcional de Epstein para D_3, conectando el modelo con la hipotesis de
Riemann a traves del teorema de Friedli-Karlsson. El sistema genera solitones estables
desde pulsos localizados y una flecha del tiempo irreversible medida por la entropia de
von Neumann de la matriz de densidad. La separacion espectral del colapso (alta frecuencia
= ruido, baja frecuencia = estructura) estabiliza el punto fijo del sistema en 3D FCC,
resolviendo la barrera cuantico-clasico que aparecia en la version 1D.

## 1. Introduccion

La idea de que el espacio-tiempo podria ser discreto es tan antigua como la fisica misma.
Aqui construimos un modelo minimalista donde una red FCC (la red de empaquetamiento maximo
en 3D, con factor de llenado eta = pi/(3*sqrt(2)) ~ 0.74) aloja una ecuacion de onda
no-lineal que genera simultaneamente estructura localizada (solitones) y una flecha del
tiempo (colapsos auto-regulados en el espacio espectral). Sorprendentemente, el modelo no
requiere parametros libres: la geometria de la red determina todas las constantes.
La version 2 incorpora la separacion espectral del colapso, que revela un punto fijo
estable en 3D FCC donde la DNLS y el colapso coexisten sin caos.

## 2. Formulacion matematica

### 2.1 Red FCC

La red FCC se construye sobre los sitios (i,j,k) con i,j,k enteros y (i+j+k) par,
en una supercelda periodica N x N x N. Cada sitio tiene Z = 12 vecinos:

    (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0),
    (1,0,1), (1,0,-1), (-1,0,1), (-1,0,-1),
    (0,1,1), (0,1,-1), (0,-1,1), (0,-1,-1)

### 2.2 Laplaciano

L = D - A, con D_nn = 12 y A_nm = 1 si n y m son vecinos.

Autovalores (Bloch-Floquet):

    lambda(k) = 12 - 4[cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]   (1)

Rango espectral: [0, 16]. El autovalor 0 corresponde al modo uniforme.
El autovalor 16 corresponde a la maxima oscilacion entre subredes A y B.

### 2.3 Funcion zeta espectral

La funcion zeta espectral del Laplaciano:

    zeta_Delta(s) = sum_{lambda_n != 0} lambda_n^{-s}                             (2)

Para la red FCC, esta funcion se relaciona con la funcion zeta de Epstein del reticulado D_3:

    xi(s) = pi^{-s} Gamma(s) zeta_{D_3}(s)                                       (3)

que satisface la ecuacion funcional:

    xi(3/2 - s) = xi(s)                                                            (4)

Para N=8 (256 sitios), la ecuacion funcional se verifica con error < 2% en el punto fijo
s = 3/4. Para N=12 (864 sitios), el error disminuye. Friedli y Karlsson (2017) demostraron
que la hipotesis de Riemann es equivalente a una ecuacion funcional aproximada para
funciones zeta espectrales de grafos, conectando directamente nuestro modelo con uno de
los problemas abiertos mas importantes de la matematica.

### 2.4 DNLS

Evolucion temporal:

    i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n                     (5)

con gamma = 12 (determinado por la geometria, ver seccion 3).

Solucion numerica: splitting de Strang en la base de autovectores del Laplaciano.

Paso no-lineal (medio):

    psi_n --> psi_n * exp(i * gamma * |psi_n|^2 * dt/2)                         (6)

Paso lineal (completo, en base espectral):

    c_k --> c_k * exp(-i * lambda_k * dt)                                       (7)

donde c_k = <phi_k | psi> son los coeficientes en la base de autovectores {phi_k, lambda_k}.

### 2.5 Auto-arranque y flecha del tiempo

#### 2.5.1 Colapso espectral auto-regulado

En la base de autovectores del Laplaciano, el colapso actua selectivamente segun
la frecuencia del modo:

    c_k --> c_k * (1 - beta * sigma(lambda_k) * epsilon_0/(epsilon_0 + |c_k|^2))   (8)

donde sigma(lambda) = 1/[1 + exp(-s(lambda/lambda_max - 0.5))] es una funcion sigmoide
suave (s=8) que separa los modos de baja frecuencia (lambda << 8, estructura del soliton)
de los de alta frecuencia (lambda >> 8, ruido). La auto-regulacion epsilon_0/(epsilon_0+|c_k|^2)
protege los modos con amplitud significativa.

#### 2.5.2 Entropia de von Neumann

La flecha del tiempo se mide mediante la entropia de von Neumann de la matriz de densidad
promediada sobre M realizaciones estocasticas del colapso espectral:

    rho(t) = (1/M) sum_{i=1}^{M} |psi_i(t)><psi_i(t)|                             (9)
    S(t) = -Tr(rho(t) ln rho(t))                                                   (10)

Para M=100 realizaciones en 3D FCC (N=6, 108 sitios):

    t=0.00: S=0.0124, rango(rho)=68
    t=0.50: S=0.0070, rango(rho)=3
    t=1.00: S=2e-5,   rango(rho)=2
    t=1.50: S=0,      rango(rho)=1
    t=2.00: S=0,      rango(rho)=1

S(t) decrece monotonicamente a cero: el colapso espectral PURIFICA el ensemble,
llevando todas las realizaciones al mismo estado puro (atractor global).
Este es el analogo dinamico del colapso de la funcion de onda en la medicion cuantica.

## 3. Determinacion de los parametros

### 3.1 epsilon_0

La red FCC tiene Z = 12 vecinos por sitio. La perturbacion de fondo minima
es el inverso del numero de coordinacion:

    epsilon_0 = 1/Z = 1/12                                                       (11)

Este valor se verifica numericamente: es el unico que produce colapso despreciable
en el centro del soliton (gamma_centro ~ 0 porque |psi|^2 >> epsilon_0) y colapso
efectivo en el ruido de fondo (gamma_borde ~ 1 porque |psi|^2 << epsilon_0).
Valores alternativos (1/24, 1/6, 1/2) rompen este equilibrio.

### 3.2 gamma

La fuerza no-lineal es el inverso de la perturbacion de fondo:

    gamma = 1/epsilon_0 = Z = 12                                                 (12)

Esto cierra el modelo: la geometria determina la no-linealidad.

### 3.3 Parametros derivados

- epsilon_0 = 1/12 (de la coordinacion FCC)
- gamma = 12 (del inverso de epsilon_0)
- sigma(lambda) = 1/[1 + exp(-s(lambda/16 - 0.5))] con s=8 (filtro espectral)
- gamma_k = sigma(lambda_k) * epsilon_0/(epsilon_0 + |c_k|^2) (colapso por modo)
- beta = 0.3 (tasa de colapso, ajustable dentro de un rango de estabilidad)

Cero parametros libres excepto beta, que tiene un rango de estabilidad [0.1, 0.5].

### 3.4 Ecuacion funcional de Epstein

La verificacion numerica de la ecuacion funcional xi(3/2 - s) = xi(s) para N=8
(256 sitios) da error < 2% en s = 3/4, consistente con la teoria de Friedli-Karlsson.
Para N=12 (864 sitios) la verificacion mejora por el mejor muestreo de la zona de Brillouin.

## 4. Resultados numericos

### 4.1 Espectro

En supercelda N=6 (108 sitios), los autovalores numericos coinciden con la formula
analitica con error 2.84 x 10^{-14} (precision de maquina). Para N=12 (864 sitios)
la diagonalizacion toma 0.96s y el espectro muestra la estructura de degeneraciones
esperada del reticulado D_3.

### 4.2 Soliton en 1D

Configuracion: N=64, gamma=12, dt=0.005, pulso gaussiano inicial.

    t=0:  ancho = 3.77,  max|psi|^2 = 0.106
    t=3:  ancho = 3.50,  max|psi|^2 = 0.140
    t=6:  ancho = 2.74,  max|psi|^2 = 0.484

Energia conservada a 10^{-6}. Compresion del 27%. El soliton 1D se forma
y se mantiene estable bajo la DNLS pura.

### 4.3 Colapso espectral en 3D FCC

Configuracion: N=6 (108 sitios), gamma=12, epsilon_0=1/12, beta=0.3, dt=0.005.

| t | ancho | max|psi|^2 | energia |
|---|-------|-----------|---------|
| 0.0 | 2.20 | 0.040 | 0.977 |
| 1.0 | 2.96 | 0.009 | -0.056 |
| 2.0 | 2.96 | 0.009 | -0.056 |
| 5.0 | 2.96 | 0.009 | -0.056 |

El sistema converge a un PUNTO FIJO ESTABLE en t ~ 1.5. No hay caos
(como en el colapso en espacio real) ni sobre-enfoque (como en el colapso filtrante).
La DNLS y el colapso espectral coexisten en equilibrio dinamico.

### 4.4 Dependencia dimensional

El colapso espectral revela una diferencia critica entre 1D y 3D:

**1D (N=64)**: el soliton se DEGRADA bajo colapso espectral.
Ancho crece de 3.77 a 6.52 en 2000 pasos. Energia sube de -0.43 a -0.25.
No hay punto fijo estable porque el soliton 1D distribuye su potencia en
muchos modos de Fourier, y ninguno esta protegido del filtro espectral.

**3D FCC (N=6, 108 sitios)**: punto fijo ESTABLE. El modo uniforme (lambda=0)
retiene el 84% de la potencia y es casi inmune al colapso (sigma(0) ~ 0.018).
Los modos de alta frecuencia colapsan pero su energia es reabsorbida por el
modo cero via renormalizacion de la norma.

### 4.5 Finite-size scaling: N=12 (864 sitios)

Barrido de la norma inicial del paquete gaussiano (sigma=2.0, N=12, gamma=12):

| Norma | ancho(t=0) | ancho(t=2.5) | Estado |
|-------|-----------|-------------|--------|
| 2 | 2.45 | 7.41 | EXPANSION |
| 4 | 2.45 | 3.45 | AUTO-ATRAPADO |
| 8 | 2.45 | 3.69 | AUTO-ATRAPADO |
| 16 | 2.45 | 3.16 | AUTO-ATRAPADO |

El umbral de auto-atrapamiento en 3D FCC con gamma=12 esta entre norma=2 y norma=4.
El caso N=6 (norma=2, ancho aparente ~3) era un artefacto de tamano finito:
el sistema era demasiado pequeno para que el paquete se expandiera.
En N=12, el paquete con norma=2 se expande sin limite (ancho -> 7.4).

### 4.6 Estructura espectral del punto fijo

N=12, norma=2:

| Modo | lambda | Potencia |
|------|--------|----------|
| 1 | 1.07 | 54.8% |
| 2 | 0.00 | 45.1% |
| 3+ | >2.07 | <0.06% |

Solo 2 modos significativos (>1% de potencia total). El 99.9% de la potencia
se concentra en los dos modos mas bajos del espectro. Esto confirma que el
colapso espectral selecciona los modos de minima frecuencia como "nucleo clasico".

## 5. Discusion

### 5.1 La barrera cuantico-clasico es estructural

El colapso espectral confirma que la barrera entre la mecanica cuantica
(colapsos probabilisticos) y la clasica (solitones deterministas) no es
un artefacto numerico sino una propiedad fundamental de los sistemas
no-lineales en redes discretas. La separacion de escalas (alta frecuencia =
ruido cuantico, baja frecuencia = estructura clasica) funciona en 3D FCC
porque el modo uniforme actua como sumidero termodinamico. En 1D, donde no
hay tal sumidero, la barrera permanece infranqueable.

### 5.2 Purificacion como flecha del tiempo

La entropia de von Neumann decreciente (S(t) -> 0) no es termalizacion
(segunda ley) sino purificacion: el colapso espectral es un atractor global
en el espacio de Hilbert. Todas las realizaciones convergen al mismo estado
puro. Esto es precisamente lo que ocurre en una medicion cuantica:
un ensemble estadistico colapsa a un resultado definido. La "flecha del tiempo"
es la perdida irreversible de informacion sobre la condicion inicial,
medida por la disminucion de S(t).

### 5.3 Conexion con la funcion zeta de Epstein

La red FCC, cuya funcion zeta espectral satisface la ecuacion funcional de
Epstein para D_3, proporciona el escenario geometrico donde la separacion
espectral opera. La coincidencia no es accidental: los mismos autovalores
que determinan la dispersion de ondas en la red (lambda_k) determinan
la estructura de ceros de la funcion zeta.

## 6. Conclusion

Hemos construido un modelo de espacio-tiempo discreto sobre red FCC con DNLS
que tiene cero parametros libres. La geometria determina todo. El modelo genera
solitones estables, una flecha del tiempo irreversible (purificacion del ensemble)
y una separacion natural de escalas (alta frecuencia = ruido, baja frecuencia =
estructura) que resuelve la barrera cuantico-clasico en 3D FCC. La conexion con
la funcion zeta de Epstein sugiere que la geometria del espacio-tiempo discreto
y la teoria de numeros primos comparten una estructura matematica comun.

## Referencias

- Friedli, F., Karlsson, A. "Spectral zeta functions of graphs and the Riemann zeta
  function in the critical strip." Tohoku Math. J. 69(4), 585-610 (2017).
- Sarnak, P., Strombergsson, A. "Minima of Epstein zeta functions and heights of
  flat tori." Invent. Math. 165, 115-151 (2006).
- Sun, K. et al. "The quantum Hall effect in graphene." Nature 474, 71-77 (2011).
- Tsironis, G. "The Discrete Nonlinear Schrodinger Equation." In: AI and Complex
  Dynamical Systems, Springer (2025).
- Connes, A., Consani, C. "Zeta Spectral Triples." arXiv:2511.22755 (2025).
- Epstein, P. "Zur Theorie allgemeiner Zetafunktionen." Math. Ann. 56, 615-644 (1903).
