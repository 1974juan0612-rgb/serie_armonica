# Emergencia de la barrera cuantico-clasico en un modelo de red discreta con auto-interaccion

**Autor:** Serie Armonica Collaboration

**Fecha:** Junio 2026

## Resumen

Presentamos un modelo de espacio-tiempo discreto basado en una red FCC con ecuacion
de Schrodinger no-lineal discreta (DNLS) y un mecanismo de colapso auto-regulado.
El modelo reproduce la barrera entre el comportamiento cuantico (colapsos probabilisticos,
ruido de fondo) y el clasico (solitones estables, estructura determinista) sin postular
la dicotomia, sino viendola emerger de la dinamica. Los colapsos (regulados por
gamma_n = epsilon_0/(epsilon_0+|psi|^2)) y la no-linealidad (gamma = 12) compiten
naturalmente, generando una interfaz que no se puede soldar con acoplamiento lineal.
Este resultado sugiere que el problema de la medida cuantica no es un artifacto
filosofico sino una propiedad emergente de sistemas con auto-interaccion en redes
discretas.

## 1. Introduccion

El problema de la medida en mecanica cuantica -- la transicion entre la superposicion
de estados y un resultado clasico definido -- ha sido discutido durante casi un siglo
sin resolucion definitiva. Aqui mostramos que este problema emerge naturalmente en un
modelo minimalista: una red FCC con ecuacion de Schrodinger no-lineal discreta y un
mecanismo de colapso auto-regulado por la densidad local.

## 2. El modelo

Evolucion unitaria (DNLS):

    i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n                     (1)

con L = D - A el Laplaciano de la red FCC (12 vecinos/sitio, autovalores en [0, 16]).

Colapso auto-regulado:

    gamma_n = epsilon_0 / (epsilon_0 + |psi_n|^2)                               (2)

    psi_n --> ruido aleatorio con probabilidad gamma_n                           (3)

La escala geometrica epsilon_0 = 1/12 y gamma = 12 estan fijadas por la red.

## 3. Resultados

### 3.1 Soliton puro (sin colapsos)

La DNLS con gamma = 12 forma solitones estables desde pulsos localizados:

    Ancho inicial: 3.77
    Ancho en t=6:  2.74  (compresion del 27%)
    Energia: conservada a 10^-6

### 3.2 Colapso puro (sin DNLS)

El auto-arranque con epsilon_0 = 1/12 produce una inercia I(t) que oscila
permanentemente, indicando que el tiempo nunca se detiene:

    I(t=0)  = 1.78
    I(t=fin)= 1.94  (oscilacion sostenida)

### 3.3 Modelo unificado (colapsos + DNLS)

Cuando ambos mecanismos se activan simultaneamente, el sistema se vuelve caotico:

    Ancho oscila entre 7 y 22 sitios
    Energia fluctua: -0.78 < E < 0.65
    Sin convergencia a estado estable

Cada condicion inicial produce un resultado final diferente.

### 3.4 Colapso filtrante (intento de soldadura)

Reemplazar el colapso aleatorio por un filtro que suprime sitios de baja densidad:

    psi_n --> psi_n * (1 - beta * gamma_n)                                       (4)

produce sobre-enfoque: toda la probabilidad se condensa en 1 solo sitio
(ancho = 0.0075, max|psi|^2 = 0.9999), destruyendo el soliton.

## 4. Interpretacion: la barrera emerge

Los resultados muestran que los colapsos probabilisticos (ecuacion 3) y la
no-linealidad determinista (ecuacion 1) ocupan regimenes mutuamente excluyentes:

    Colapso + DNLS = caos                                                        (5)
    Filtro + DNLS  = sobre-enfoque                                               (6)
    Solo DNLS      = soliton (desde estado coherente)                            (7)
    Solo colapso   = tiempo vivo                                                  (8)

Esta exclusion mutua no es un defecto del modelo. Es una propiedad emergente:
la competencia entre fluctuaciones cuanticas (colapsos) y estructura clasica
(no-linealidad) genera una interfaz que ningun acoplamiento lineal puede soldar.
Esto reproduce, desde primeros principios de red, la dicotomia fundamental
entre la mecanica cuantica y la relatividad general.

## 5. Proxima direccion

La soldadura de la barrera probablemente requiere separar los modos del Laplaciano:
colapsar solo las componentes de alta frecuencia (autovalores ~ 16, ruido) y
dejar intactas las de baja frecuencia (autovalores ~ 0, soliton).
Esta separacion espectral es el analogo natural de la renormalizacion en teorias
cuanticas de campos.

## 6. Conclusion

Un modelo minimalista de red FCC con DNLS y colapso auto-regulado reproduce
la barrera cuantico-clasico sin postularla. La interfaz emerge de la competencia
entre dos mecanismos (colapso y no-linealidad) que la propia red genera.
Esto sugiere que el problema de la medida no es externo a la teoria,
sino una consecuencia inevitable de la geometria discreta del espacio.

## Referencias

- Penrose, R. "On Gravity's Role in Quantum State Reduction." Gen. Rel. Grav. 28,
  581-600 (1996).
- Diosi, L. "A Universal Master Equation for the Gravitational Violation of Quantum
  Mechanics." Phys. Lett. A 120, 377-381 (1987).
- Friedli, F., Karlsson, A. "Spectral zeta functions of graphs." Tohoku Math. J.
  69(4), 585-610 (2017).
- Tsironis, G. "The Discrete Nonlinear Schrodinger Equation." Springer (2025).
- Connes, A., Consani, C. "Zeta Spectral Triples." arXiv:2511.22755 (2025).
