# Horizontes de Investigación: Tres Preguntas Abiertas

**Modelo:** FCC-DNLS con colapso espectral y cero parámetros libres
**Autor:** Serie Armónica Collaboration
**Base:** Manuscrito v3 (β = (1-η) + ε₀/2, mapa de Kraus causal)
**Fecha:** Junio 2026

---

## Introducción

Todo marco teórico de frontera se valida por su poder de generar
preguntas que antes no podían formularse. El modelo FCC-DNLS
resuelve problemas históricos (colapso cuántico sin parámetros
libres, dimensionalidad del espacio, origen de la masa como
auto-atrapamiento), pero también abre tres horizontes
radicalmente nuevos. A continuación formalizamos cada uno con
predicciones comprobables y programas de investigación.

---

## A. Pregunta Cosmológica: ¿Es el universo temprano un efecto
   de la auto-interferencia espectral de un gran solitón?

### El problema

El telescopio James Webb (JWST) detecta galaxias masivas
($M_\star \sim 10^{10-11} M_\odot$) y agujeros negros
supermasivos ($M_{\rm BH} \sim 10^7 M_\odot$) a
redshifts $z > 10$, cuando el universo tenía menos de
500 millones de años. En el modelo $\Lambda$CDM, esto requiere
eficiencias de formación estelar imposibles ($\epsilon > 1$)
o semillas de agujeros negros primordiales.

### La respuesta del modelo

En el formalismo FCC-DNLS, el universo observable sería un
solitón espectral en una red FCC infinita. La conexión entre
distancia cosmológica y desplazamiento al rojo no es la
expansión métrica FLRW, sino el **gradiente espectral del
perfil del solitón**:

$$\frac{\delta\lambda}{\lambda} \approx \frac{r^2}{2\xi^2}$$

donde $r$ es la distancia al centro del solitón (nuestra
posición observacional) y $\xi \approx 2.96$ es el ancho
característico.

JWST ve galaxias masivas a $z > 10$ no porque existieran
temprano en el tiempo, sino porque el **mapeo distancia-redshift
no es lineal**. Un objeto al borde del solitón ($r \sim \xi$)
sufre un corrimiento espectral $\delta\lambda/\lambda \sim 0.5$,
que en el marco FLRW se interpreta como $z \sim 1$ pero que
en realidad corresponde a una distancia física mucho menor.

### Predicciones falsables

| # | Predicción | Observable | Falsación |
|---|------------|------------|-----------|
| A1 | La relación distancia-redshift se desvía de FLRW para $z > 6$ siguiendo una ley cuadrática, no lineal | Curva de Hubble $H(z)$ de SNIa + JWST | Si $H(z)$ sigue FLRW exactamente hasta $z=15$ |
| A2 | Las galaxias a $z > 10$ tienen morfologías "espectralmente replicadas": misma estructura interna que galaxias a $z \sim 2$ pero proyectadas con distorsión radial | Morfología de galaxias JWST vs. HST CANDELS | Si las galaxias $z>10$ muestran morfologías únicas sin análogo local |
| A3 | El fondo cósmico de microondas (CMB) tiene una modulación dipolar residual del perfil del solitón, no atribuible al movimiento del sistema solar | Dipolo CMB en mapas de Planck | Si el dipolo CMB se explica completamente por $v_{\rm solar} = 370$ km/s |

### Programa de investigación

1. Derivar el mapeo explícito entre la métrica inducida por el
   solitón FCC-DNLS y la métrica FLRW: $g_{\mu\nu}^{\rm soliton}$
   como función del perfil $|\psi(r)|^2$ y el ancho $\xi$.
2. Simular catálogos de galaxias sintéticos en el espacio de
   fase FCC y comparar con JWST, Euclid y Roman.
3. Predecir la firma del "eco espectral": objetos replicados a
   intervalos periódicos en $z$ correspondientes a la
   periodicidad de la red $L = N a$.

---

## B. Pregunta Cuántica: ¿Es la no-localidad cuántica una
   ilusión de la base espacial?

### El problema

Las violaciones de la desigualdad de Bell (Aspect 1982,
Hensen 2015) demuestran que las correlaciones cuánticas no
pueden explicarse por variables ocultas locales. Esto se
interpreta como "acción fantasmal a distancia" (Einstein) o
como no-localidad irreducible.

### La respuesta del modelo

En el modelo FCC-DNLS, el colapso es **local en la base
espectral** (autovectores del Laplaciano). La transformación
entre la base espectral y la base real (sitios de la red)
es la transformada de Fourier discreta:

$$\psi(n) = \sum_k c_k \phi_k(n), \quad \phi_k(n) \sim e^{i k \cdot n}$$

Un operador multiplicativo en $k$-espacio es una convolución
en $x$-espacio:

$$\psi'(n) = \psi(n) - \beta \sum_m K(n-m) \psi(m)$$

$$K(\Delta n) = \sum_k \sigma(\lambda_k) \frac{\varepsilon_0}{\varepsilon_0+|c_k|^2} e^{i k \cdot \Delta n}$$

El núcleo $K(\Delta n)$ tiene **soporte compacto** determinado
por el filtro sigmoide: $K(\Delta n) \approx 0$ para
$|\Delta n| > \Lambda \approx 8/\Delta k \sim 2$ sitios.

La "no-localidad" cuántica emerge porque:

1. Las mediciones se realizan en base real (detectores en
   posiciones específicas)
2. El colapso es local en base espectral
3. La transformada de Fourier proyecta esta localidad espectral
   como correlación instantánea en espacio real

Esto NO viola causalidad porque la dinámica del colapso es
Markoviana (mapa de Kraus, sin memoria) y la tasa de
propagación en espacio real satisface $v_{\max} \leq
\Delta x / \Delta t \approx 0.5$ sitios/unidad, muy por debajo
del cono de luz.

Las correlaciones de Bell son entonces **correlaciones
espectrales pre-existentes** en la función de onda conjunta
del sistema-lattice, que el colapso revela pero no transmite.

### Predicciones falsables

| # | Predicción | Observable | Falsación |
|---|------------|------------|-----------|
| B1 | La violación de Bell decae con la distancia como $K(\Delta n)$ — correlaciones no nulas sólo hasta $|\Delta n| \lesssim 2\xi$ | Experimentos Bell con separación variable en redes ópticas FCC | Si la violación persiste sin decaimiento a distancias $>4\xi$ |
| B2 | El colapso espectral induce una fase geométrica medible en interferometría, proporcional a $\beta\varepsilon_0$ | Interferometría atómica en redes ópticas | Si la fase observada es nula |
| B3 | El teorema de Bell no se viola en la base espectral: existe una descripción local oculta en $k$-espacio | Tomografía de procesos cuánticos en base de Fourier | Si se demuestra que NO existe modelo local en $k$-espacio |

### Programa de investigación

1. Derivar la cota de CHSH en función de la distancia entre
   sitios en la red FCC, usando el núcleo $K(\Delta n)$.
2. Construir explícitamente la variable oculta local en
   $k$-espacio: $\lambda_k = \{\arg(c_k), |c_k|^2\}$ y mostrar
   que la estadística de Bell es reproducible.
3. Proponer un experimento de Bell en redes ópticas FCC
   (similar a los de Monz et al., 2009, pero con colapso
   controlado por la profundidad del potencial óptico).

---

## C. Pregunta Matemática: ¿Cuál es el papel exacto de los
   números primos en el espacio-tiempo discreto?

### El problema

La función zeta de Riemann $\zeta(s) = \sum_n n^{-s}$ codifica
la distribución de números primos. La hipótesis de Riemann
(los ceros no triviales están en $\operatorname{Re}(s)=1/2$)
sigue sin demostración. La función zeta de Epstein para la red
$D_3$ satisface $\xi(3/2-s) = \xi(s)$, con línea crítica en
$\operatorname{Re}(s) = 3/4$.

### La respuesta del modelo

El Laplaciano de la red FCC genera autovalores $\lambda_k$
cuya distribución espectral determina:

1. La dispersión de ondas (propagación de solitones)
2. La función zeta espectral $\zeta_\Delta(s) = \sum_{\lambda_k \neq 0} \lambda_k^{-s}$
3. La densidad de estados $\rho(\lambda) = \frac{1}{N} \sum_k \delta(\lambda - \lambda_k)$
4. El umbral de auto-atrapamiento (norma $= 4$)
5. La estabilidad del atractor global (33 modos de baja frecuencia)

Conexión profunda: los autovalores del Laplaciano FCC son

$$\lambda(k) = 12 - 4[\cos k_x\cos k_y + \cos k_x\cos k_z + \cos k_y\cos k_z]$$

Para $k$ pequeños, $\lambda \sim |k|^2$. La función zeta
espectral $\zeta_\Delta(s)$ y la función zeta de Epstein
$\zeta_{D_3}(s)$ están relacionadas por el teorema de
Friedli-Karlsson:

$$\xi_\Delta(s) = \xi_{D_3}(s) + \text{términos de red finita},$$

$$\xi(s) = \pi^{-s}\Gamma(s)\zeta_{D_3}(s),\quad \xi(3/2-s) = \xi(s).$$

El error numérico $<2\%$ en $s=3/4$ para $N=8$ (256 sitios)
confirma que la conexión es estructural.

### La hipótesis de trabajo

La distribución de autovalores $\lambda_k$ del Laplaciano FCC
codifica información aritmética sobre las representaciones de
enteros como suma de tres cuadrados (el número de
representaciones está dado por los coeficientes de la función
theta de Jacobi $\vartheta_3(q)^3$). Esta información es
**legible por el colapso espectral**: el filtro
$\sigma(\lambda)$ selecciona modos con $\lambda \ll 8$, y la
estructura fina de $\rho(\lambda)$ cerca de $\lambda=0$ está
determinada por la aritmética del retículo $D_3$.

Si el umbral de auto-atrapamiento (norma $=4$) corresponde a
un autovalor $\lambda_4$ con propiedades aritméticas especiales
(por ejemplo, $\lambda_4$ siendo un número altamente compuesto
en el espectro), entonces la **existencia de materia clásica
estable sería contingente a la distribución de números primos**:
una forma radical de la hipótesis de Riemann aplicada a la
física.

### Predicciones falsables

| # | Predicción | Observable | Falsación |
|---|------------|------------|-----------|
| C1 | El error en $\xi(3/2-s) = \xi(s)$ escala como $O(N^{-\alpha})$ con $\alpha \approx 1/2$ | Medición numérica para $N=4,6,8,10,12,16$ | Si $\alpha < 0.1$ (no converge) |
| C2 | El umbral de auto-atrapamiento (norma $=4$) es un autovalor aritméticamente distinguido del Laplaciano FCC | Factorización del polinomio característico de $L$ | Si norma $4$ no corresponde a ningún autovalor especial |
| C3 | Los 33 canales de baja frecuencia se organizan en multipletes según el grupo de Weyl de $D_3$, con degeneraciones dadas por coeficientes de representaciones enteras | Descomposición del espectro bajo el grupo de simetrías FCC | Si las degeneraciones no corresponden a caracteres de $D_3$ |

### Programa de investigación

1. Demostrar que el teorema de Friedli-Karlsson se extiende
   a la red finita con error controlado: $|\xi_\Delta(s) - \xi_{D_3}(s)| \leq CN^{-1/2}$.
2. Calcular la distribución de autovalores cerca de $\lambda=0$
   y mostrar que $\rho(\lambda) \sim \lambda^{D/2-1} \times$
   (factores aritméticos) para $D=3$.
3. Investigar si la equivalencia espectral del colapso
   (selección de modos $\lambda \approx 0$) es un análogo
   físico del "cernido" de números primos en la criba de
   Eratóstenes, con los modos colapsados correspondiendo a
   números compuestos y los modos protegidos a números primos.

---

## Conclusión: Tres Programas de Investigación

```
PR: Cosmología (JWST)    PR: No-localidad (Bell)    PR: Números primos (zeta)
       |                         |                          |
       v                         v                          v
   Mapeo solitón-FLRW     Límite CHSH vs distancia    Error zeta ~ N^{-1/2}
   Catálogos sintéticos   Variable oculta en k-space   Umbral de norma = 4
   Eco espectral          Fase geométrica βε₀         Multipletes D₃
```

Cada programa genera predicciones falsables explícitas con
tecnología actual (redes ópticas FCC, JWST, interferometría
atómica). El modelo FCC-DNLS no es una especulación
filosófica: es un marco computacional y matemático que puede
ser verificado o refutado en un horizonte de 3-5 años.

---
*Serie Armónica Collaboration*
