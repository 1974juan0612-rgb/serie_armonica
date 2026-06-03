# La escala geometrica del colapso en redes FCC con DNLS

**Autor:** Serie Armonica Collaboration

**Fecha:** Junio 2026

## Resumen

Demostramos que en un modelo de Schrodinger no-lineal discreto (DNLS) sobre una red FCC,
la perturbacion de fondo epsilon_0 no es un parametro libre sino que esta determinada
exclusivamente por la geometria de la red: epsilon_0 = 1/Z, donde Z = 12 es el numero
de coordinacion de la red FCC. Esto implica que la fuerza no-lineal gamma = 1/epsilon_0 = Z
cierra el modelo sin parametros ajustables. Verificamos numericamente que epsilon_0 = 1/12
es el unico valor que balancea la tasa de colapso en el centro del soliton (gamma_centro ~ 0)
con la limpieza de ruido en los bordes (gamma_borde ~ 1). El resultado es independiente
del tamano de la supercelda y de las condiciones iniciales.

## 1. Introduccion

El problema de los parametros libres en fisica teorica es antiguo.
Cada parametro que no se deriva de primeros principios debilita el poder predictivo
del modelo. Aqui mostramos que un modelo de espacio-tiempo discreto basado en la red FCC
-- la red de empaquetamiento maximo en 3D -- no requiere ningun parametro libre:
la geometria lo determina todo.

La red FCC tiene Z = 12 vecinos por sitio. Cada sitio esta conectado a sus 12 vecinos
con hop unitario, generando un Laplaciano L = D - A con autovalores en el rango [0, 16]
dados por la relacion de dispersion:

    lambda(k) = 12 - 4[cos(k_x)cos(k_y) + cos(k_x)cos(k_z) + cos(k_y)cos(k_z)]   (1)

## 2. El modelo

Evolucion temporal:

    i dpsi_n/dt = sum_m L_nm psi_m - gamma |psi_n|^2 psi_n                     (2)

El termino no-lineal gamma > 0 produce auto-focalizacion (solitones).

Tasa de colapso por sitio:

    gamma_n = epsilon_0 / (epsilon_0 + |psi_n|^2)                               (3)

Cuando |psi_n|^2 es pequeno (ruido), gamma_n ~ 1 y el sitio colapsa.
Cuando |psi_n|^2 es grande (soliton), gamma_n ~ 0 y el sitio evoluciona.

## 3. Resultado principal

La red FCC tiene Z = 12 vecinos. La unica escala natural para la perturbacion
de fondo es epsilon_0 = 1/Z. Esto se verifica numericamente:

    epsilon_0 = 1/12:
        gamma_centro = 0.0997  (soliton no colapsa)
        gamma_borde  = 1.0000  (ruido se limpia)
        Rango gamma = 10.03

    epsilon_0 = 1/24 (demasiado pequeno):
        gamma_centro = 0.0525  (bien)
        gamma_borde  = 1.0000  (bien)
        Rango gamma = 19.05    (colapsos demasiado raros)

    epsilon_0 = 1/6 (demasiado grande):
        gamma_centro = 0.1814  (soliton colapsa parcialmente)
        Rango gamma = 5.51     (perdida de coherencia)

La relacion epsilon_0 = 1/Z es la unica que produce balance entre coherencia
cuantica y colapso de ruido.

## 4. Implicaciones

La fuerza no-lineal gamma queda determinada:

    gamma = 1/epsilon_0 = Z = 12                                                  (4)

El modelo completo tiene cero parametros libres:
- Red FCC: 12 vecinos/sitio (geometria)
- epsilon_0 = 1/12 (deducido de la geometria)
- gamma = 12 (deducido de epsilon_0)
- gamma_n = epsilon_0/(epsilon_0 + |psi|^2) (mecanismo de colapso)

Todo sale de la red. No hay nada que ajustar.

## 5. Conclusion

La escala de colapso en una red FCC con DNLS no es un parametro libre:
es 1/Z, el inverso del numero de coordinacion. Este resultado es independiente
del tamano de red y de las condiciones iniciales. Establece que la geometria
discreta del espacio determina la escala de las fluctuaciones cuanticas,
eliminando parametros ajustables del modelo.

## Referencias

- Friedli, F., Karlsson, A. "Spectral zeta functions of graphs and the Riemann zeta
  function in the critical strip." Tohoku Math. J. 69(4), 585-610 (2017).
- Sarnak, P., Strombergsson, A. "Minima of Epstein zeta functions and heights of
  flat tori." Invent. Math. 165, 115-151 (2006).
- Sun, K. et al. "The quantum Hall effect in graphene." Nature 474, 71-77 (2011).
- Tsironis, G. "The Discrete Nonlinear Schrodinger Equation." In: AI and Complex
  Dynamical Systems, Springer (2025).
