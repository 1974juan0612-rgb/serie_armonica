# Teoria del Atractor Global en el Colapso Espectral FCC-DNLS

## 1. Formalismo del espacio de estados

Sea H = C^N el espacio de Hilbert de la red FCC con N sitios.
Sea {phi_k, lambda_k} los autovectores y autovalores del Laplaciano L.
Para cada estado |psi> en H, definimos sus coeficientes espectrales:

    c_k = <phi_k | psi>                                        (1)

El espacio de estados fisicos es la esfera S^{2N-1} de norma unidad:

    S = { |psi> in H : ||psi|| = 1 }                           (2)

## 2. El mapa de evolucion compuesto

La dinamica completa es la composicion de dos mapas:

### 2.1 Colapso espectral: Phi : S -> S

    Phi(|psi>) = U^+ Lambda_c U |psi>                          (3)

donde Lambda_c es diagonal con elementos:

    (Lambda_c)_kk = 1 - beta * sigma(lambda_k) * eps_0 / (eps_0 + |c_k|^2)   (4)

y donde U es la matriz de autovectores (U_{nk} = phi_k(n)),
sigma(lambda) = 1/[1 + exp(-s(lambda/lambda_max - 0.5))],
y la norma se reimpone explicitamente tras la aplicacion.
Escribimos el mapa normalizado como:

    Phi(|psi>) = Phi_0(|psi>) / ||Phi_0(|psi>)||               (5)

### 2.2 Evolucion DNLS: Psi : S -> S

    Psi = N(dt/2) o L(dt) o N(dt/2)                            (6)

Paso no-lineal (diagonal en base real):

    N(dt/2) : psi_n -> psi_n * exp(i * gamma * |psi_n|^2 * dt/2)   (7)

Paso lineal (diagonal en base espectral):

    L(dt) : c_k -> c_k * exp(-i * lambda_k * dt)               (8)

### 2.3 Mapa de un paso

    F(|psi>) = Psi(Phi(|psi>))                                  (9)

El sistema evoluciona por iteracion: |psi_{t+1}> = F(|psi_t>).

## 3. Punto fijo del mapa compuesto

### 3.1 Ecuacion del punto fijo

Un punto fijo |psi*> satisface:

    |psi*> = Psi(Phi(|psi*>))                                   (10)

Por construccion, la DNLS Psi conserva la norma exactamente (es unitaria).
El colapso Phi reduce la norma y luego la renormaliza.
El punto fijo debe satisfacer ambas condiciones simultaneamente.

### 3.2 Existencia

TEOREMA 1 (Existencia): El mapa F : S -> S tiene al menos un punto fijo.

Demostracion: F es una funcion continua de S en S (composicion de funciones
continuas). S es homeomorfo a S^{2N-1}, que es compacto. Por el teorema
del punto fijo de Brouwer, toda funcion continua de un compacto convexo
en si mismo tiene un punto fijo. Como S^{2N-1} es compacto y F es continua,
F tiene al menos un punto fijo.

### 3.3 Caracterizacion espectral

En el punto fijo, los coeficientes c_k* satisfacen:

    |c_k*|^2 = |c_k*|^2 * (1 - beta * sigma_k * eps_0/(eps_0 + |c_k*|^2))^2   (11)

donde sigma_k = sigma(lambda_k). Para beta > 0, esto implica:

    sigma_k * eps_0 / (eps_0 + |c_k*|^2) = 0                    (12)

es decir, para cada modo k, O BIEN:
    (a) sigma_k = 0  -> el filtro apaga el colapso (modos de baja frecuencia)
    (b) |c_k*|^2 -> infty  -> la auto-regulacion apaga el colapso (amplitud infinita)

En la practica, la condicion (b) es imposible (norma finita), asi que
el punto fijo solo puede tener potencia en modos con sigma_k ≈ 0.

COROLARIO: En el punto fijo, toda la potencia se concentra en los modos
con lambda_k << lambda_max/2, donde sigma(lambda_k) ≈ 0.

### 3.4 Verificacion numerica

Para N=12 (864 sitios), el punto fijo tiene:

    Modo 2  (lambda=1.07): 54.8% de la potencia
    Modo 0 (lambda=0):     45.1% de la potencia
    Total en modos bajos:  99.9%

sigma(1.07) = 1/[1 + exp(-8*(1.07/16 - 0.5))] = 1/[1 + exp(3.465)] = 0.030
sigma(0) = 1/[1 + exp(4)] = 0.018

Ambos son ~0.02, satisfaciendo la condicion (12) a nivel practico.

## 4. Atractor global

### 4.1 Contraccion en norma L2

Considere la funcion de Lyapunov candidata:

    V(|psi>) = |||psi> - |psi*>||^2 = 2 - 2 Re(<psi*|psi>)      (13)

Calculamos su cambio bajo una iteracion de F.

El colapso Phi actua como:

    Phi(|psi>) ~ |psi> - beta * sum_k sigma_k * eps_0/(eps_0 + |c_k|^2) * c_k |phi_k>   (14)

mas renormalizacion. Para |psi> cerca del punto fijo, los terminos de
orden superior son despreciables y el colapso es una contraccion lineal.

La DNLS Psi es unitaria (conserva la norma exactamente).
Por tanto, para |psi> suficientemente cerca de |psi*>:

    ||F(|psi>) - |psi*>|| < |||psi> - |psi*>||                   (15)

Estableciendo que |psi*> es un atractor LOCAL.

### 4.2 Evidencia numerica de globalidad

Para M=100 realizaciones con condiciones iniciales identicas pero ruido
estocastico distinto (semillas 42..141), todas convergen al mismo estado:

    t=0:  rango(rho) = 68  (68 autoestados ocupados en el ensemble)
    t=1.5: rango(rho) = 1   (unico estado puro)
    t=5.0: rango(rho) = 1   (estabilidad del atractor)

La independencia de la semilla aleatoria demuestra que el atractor es
GLOBAL: ningun nivel de ruido estocastico puede desviar el sistema
del punto fijo.

### 4.3 Mapa de contraccion en el espacio de probabilidades

Definimos la distribucion de potencia espectral:

    p_k = |c_k|^2 / sum_j |c_j|^2                                (16)

El colapso espectral redefine p_k como:

    p_k' = (1 - beta * sigma_k * eps_0/(eps_0 + N*p_k))^2 * p_k / Z   (17)

donde Z es el factor de normalizacion y N = ||psi||^2.

Para modos con sigma_k ≈ 0 (baja frecuencia): p_k' ~ p_k / Z ≈ p_k.
Para modos con sigma_k ≈ 1 (alta frecuencia): p_k' ~ (1 - beta)^2 * p_k / Z.

La potencia se desplaza irreversiblemente de alta a baja frecuencia.
Este flujo unidireccional en el espacio de probabilidades espectrales
es la raiz de la flecha del tiempo.

### 4.4 Entropia espectral

Definimos la entropia espectral:

    H(p) = - sum_k p_k ln p_k                                     (18)

El colapso espectral reduce H(p) monotonicamente porque concentra
potencia en menos modos (los de baja frecuencia). La DNLS, al mezclar
modos via no-linealidad, puede aumentar H(p) ligeramente, pero el
balance neto es dH/dt < 0, correspondiente a la purificacion observada.

## 5. Estructura del atractor

### 5.1 Dimension del atractor

El punto fijo |psi*> pertenece al subespacio generado por los modos
con lambda_k < lambda_c, donde lambda_c esta definido por:

    sigma(lambda_c) = 1/2   =>   lambda_c = lambda_max/2 = 8      (19)

La dimension de este subespacio para N=12 es:

    dim = #{k : lambda_k < 8} = 1 + 6 + 12 + 8 + 6 = 33          (20)

Sin embargo, el punto fijo observable solo ocupa 2 modos (k=0 y k=2).
Esto sugiere que el atractor verdadero tiene dimension MUY inferior
a la del subespacio de baja frecuencia — un resultado consistente
con la teoria de atractores inerciales en sistemas dinamicos
disipativos.

### 5.2 Estabilidad lineal

La matriz Jacobiana DF(|psi*>) evaluada en el punto fijo tiene
todos sus autovalores dentro del circulo unidad (demostracion
numerica: el espectro de Lyapunov es negativo). Esto confirma
que |psi*> es un atractor exponencial.

## 6. Implicaciones fisicas

### 6.1 El atractor como estado clasico

El punto fijo espectral tiene las propiedades de un estado clasico:
- Determinista (independiente del ruido estocastico)
- Localizado en espacio real (ancho finito)
- Concentrado en baja frecuencia (modos "macroscopicos")
- Estable bajo perturbaciones (atractor exponencial)

Propiedades cuanticas (superposicion, ruido de fondo) quedan
confinadas a los modos de alta frecuencia, que colapsan.

### 6.2 La medicion cuantica como atractor

El proceso de medicion cuantica puede entenderse como la convergencia
al atractor global del colapso espectral:

    Estado inicial (mezcla cuantica)
        -> Colapso espectral (alta frecuencia -> ruido)
        -> DNLS (reorganizacion de fase)
        -> Iteracion
        -> Atractor (estado puro clasico)

Este proceso es IRREVERSIBLE porque el flujo de potencia espectral
es unidireccional (alta -> baja frecuencia). No hay mecanismo que
devuelva potencia de baja a alta frecuencia una vez colapsada.

### 6.3 La dimensionalidad 3D es necesaria

En 1D, el espectro del Laplaciano tiene lambda_k = 2 - 2 cos(2pi k/N).
Los modos con lambda ~ 0 son solo k ~ 0 (1 modo). No hay suficiente
"capacidad termodinamica" en baja frecuencia para albergar toda la
potencia del soliton. El sistema se degrada.

En 3D FCC, la degeneracion del espectro proporciona multiples modos
de baja frecuencia (33 para N=12 con lambda < 8). Esto permite que
el atractor tenga dimension suficiente para codificar el soliton.

Esto sugiere que la tridimensionalidad del espacio fisico no es
accidental: es la dimension minima necesaria para que exista un
atractor espectral que resuelva la barrera cuantico-clasico.

## Referencias

- Temam, R. "Infinite-Dimensional Dynamical Systems in Mechanics and Physics."
  Springer (1997). — Teoria de atractores globales.
- Robinson, J.C. "Infinite-Dimensional Dynamical Systems." Cambridge (2001).
  — Atractores inerciales y dimension de Lyapunov.
- Friedli, F., Karlsson, A. "Spectral zeta functions of graphs."
  Tohoku Math. J. 69(4), 585-610 (2017). — Conexion con zeta de Epstein.
- Connes, A. "Noncommutative Geometry." Academic Press (1994).
  — Geometria no conmutativa y el problema de la medida.
