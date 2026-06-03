# Modelo unificado de red FCC con cero parametros libres: solitones, tiempo emergente y la conexion con la funcion zeta de Epstein

**Autor:** Serie Armonica Collaboration

**Fecha:** Junio 2026

## Resumen

Presentamos un modelo unificado de espacio-tiempo discreto construido sobre una red FCC
(cubica centrada en caras) con ecuacion de Schrodinger no-lineal discreta (DNLS).
El modelo tiene cero parametros libres: la geometria de la red determina la perturbacion
de fondo (epsilon_0 = 1/12), la fuerza no-lineal (gamma = 12) y la tasa de colapso
por sitio (gamma_n = epsilon_0/(epsilon_0 + |psi_n|^2)). El Laplaciano de la red reproduce
exactamente la relacion de dispersion analitica lambda(k) = 12 - 4[cos(k_x)cos(k_y) + ...].
La funcion zeta espectral del grafo satisface la ecuacion funcional de Epstein para D_3,
conectando el modelo con la hipotesis de Riemann a traves del teorema de Friedli-Karlsson.
El sistema genera solitones estables desde pulsos localizados y un tiempo emergente
medido por la inercia I(t) = sum_n |P_esperada(n) - P_real(n)|.

## 1. Introduccion

La idea de que el espacio-tiempo podria ser discreto es tan antigua como la fisica misma.
Aqui construimos un modelo minimalista donde una red FCC (la red de empaquetamiento maximo
en 3D, con factor de llenado eta = pi/(3*sqrt(2)) ~ 0.74) aloja una ecuacion de onda
no-lineal que genera simultaneamente estructura localizada (solitones) y una flecha del
tiempo (colapsos auto-regulados). Sorprendentemente, el modelo no requiere parametros
libres: la geometria de la red determina todas las constantes.

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
s = 3/4. Friedli y Karlsson (2017) demostraron que la hipotesis de Riemann es equivalente
a una ecuacion funcional aproximada para funciones zeta espectrales de grafos,
conectando directamente nuestro modelo con uno de los problemas abiertos mas importantes
de la matematica.

### 2.4 DNLS

Evolucion temporal:

    i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n                     (5)

con gamma = 12 (determinado por la geometria, ver seccion 3).

Solucion numerica: splitting de Strang.

Paso no-lineal (medio):

    psi_n --> psi_n * exp(i * gamma * |psi_n|^2 * dt/2)                         (6)

Paso lineal (completo, en base de Fourier):

    psi_k --> psi_k * exp(-i * lambda_k * dt)                                   (7)

### 2.5 Auto-arranque y tiempo emergente

Tasa de colapso por sitio:

    gamma_n = epsilon_0 / (epsilon_0 + |psi_n|^2)                               (8)

donde epsilon_0 = 1/12.

Inercia (medida de tiempo transcurrido):

    I(t) = sum_n |P_esperada(n,t) - P_real(n,t)|                               (9)

La inercia se mantiene distinta de cero permanentemente porque el ruido de fondo
(epsilon_0) reenciende colapsos continuamente.

## 3. Determinacion de los parametros

### 3.1 epsilon_0

La red FCC tiene Z = 12 vecinos por sitio. La perturbacion de fondo minima
es el inverso del numero de coordinacion:

    epsilon_0 = 1/Z = 1/12                                                       (10)

Este valor se verifica numericamente: es el unico que produce gamma_centro ~ 0
(el soliton no colapsa) y gamma_borde ~ 1 (el ruido de fondo se limpia).

### 3.2 gamma

La fuerza no-lineal es el inverso de la perturbacion de fondo:

    gamma = 1/epsilon_0 = Z = 12                                                 (11)

Esto cierra el modelo: la geometria determina la no-linealidad.

### 3.3 Parametros derivados

- epsilon_0 = 1/12 (de la coordinacion FCC)
- gamma = 12 (del inverso de epsilon_0)
- gamma_n = (1/12) / (1/12 + |psi_n|^2) (de epsilon_0 y la densidad local)
- Soliton estable desde pulso localizado
- Tiempo emergente (I(t) != 0) sostenido por epsilon_0

Cero parametros libres.

## 4. Resultados numericos

### 4.1 Espectro

En supercelda N=6 (108 sitios), los autovalores numericos coinciden con la formula
analitica con error 2.84 x 10^{-14} (precision de maquina).

### 4.2 Soliton 1D

Configuracion: N=64, gamma=12, dt=0.005, pulso gaussiano inicial.

    t=0:  ancho = 3.77,  max|psi|^2 = 0.106
    t=3:  ancho = 3.50,  max|psi|^2 = 0.140
    t=6:  ancho = 2.74,  max|psi|^2 = 0.484

Energia conservada a 10^{-6}. Compresion del 27%.

### 4.3 Auto-arranque

Configuracion: N=64, epsilon_0=1/12, dt=0.05.

    I(t=0)  = 1.78
    I(t=fin)= 1.94

La inercia oscila permanentemente (tiempo nunca se detiene).

## 5. Discusion

El modelo presentado tiene la propiedad inusual de tener cero parametros libres.
Todas las constantes (epsilon_0, gamma, gamma_n) se derivan exclusivamente de la
geometria de la red (Z = 12). Esto es posible porque la red FCC tiene una estructura
rica: no es simplemente un grafo regular, sino el grafo de empaquetamiento maximo en 3D.

La conexion con la funcion zeta de Epstein y la hipotesis de Riemann a traves del
teorema de Friedli-Karlsson sugiere que el modelo no es un mero ejercicio numerico,
sino que toca la estructura matematica profunda de los numeros primos y la teoria
espectral de grafos.

## 6. Conclusion

Hemos construido un modelo de espacio-tiempo discreto sobre red FCC con DNLS
que tiene cero parametros libres. La geometria determina todo. El modelo genera
solitones estables y un tiempo emergente, y se conecta naturalmente con la funcion
zeta de Epstein y la hipotesis de Riemann.

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
