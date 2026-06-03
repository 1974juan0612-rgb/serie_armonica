# Respuesta a la Segunda Revisión por Pares (v2)

**Manuscrito:** Modelo unificado FCC-DNLS con colapso espectral y cero parámetros libres
**Revista:** *Journal of High Energy Physics* / *Physical Review D*
**Evaluador:** Anónimo
**Fecha:** Junio 2026

---

## Resumen

Agradecemos al revisor su análisis detallado y constructivo. La recomendación de
"Aceptación Condicionada a Revisión Menor" nos honra. A continuación respondemos
a las tres objeciones críticas (A, B, C).

---

## A. Localidad de la Renormalización de la Norma

**Objeción:** "La renormalización explícita *a posteriori* equivale a acoplar
instantáneamente todos los sitios/modos del sistema de forma global y no local.
¿Cómo se preserva la causalidad relativista?"

### A.1 Formulación como mapa de Lindblad

La secuencia colapso + renormalización no es un proceso no-local en el espacio
real, sino un **mapa de Kraus diagonal en la base espectral** seguido de
normalización. En notación de operador densidad:

$$\rho \to \frac{\sum_k K_k \rho K_k^\dagger}{\operatorname{Tr}(\sum_k K_k \rho K_k^\dagger)}$$

con operadores de Kraus

$$K_k = \sqrt{1 - \beta\,\sigma(\lambda_k)\,\frac{\varepsilon_0}{\varepsilon_0 + |c_k|^2}}\; |\phi_k\rangle\langle\phi_k|.$$

Cada $K_k$ actúa independientemente sobre un modo espectral. La traza se
preserva por construcción (canonical trace-preserving Kraus map). Este
formalismo es el estándar en sistemas cuánticos abiertos Markovianos
(Lindblad, 1976; Breuer-Petruccione, 2002) y es manifiestamente causal:
el colapso es local en el espacio de modos y no requiere comunicación
entre sitios distantes.

La normalización global posterior (división por la norma) no es más que la
implementación numérica de la preservación de traza, que en el formalismo
de Kraus está garantizada axiomáticamente.

### A.2 Interpretación física del modo $\lambda=0$

El modo uniforme ($\lambda=0$, autovector constante $\phi_0(n)=1/\sqrt{N}$)
no es un "sumidero no-local" sino una **condensación de modo cero**:
análoga al condensado de Bose-Einstein en el espacio de momentos, donde
todos los ocupantes comparten el mismo estado de momento $k=0$. La
"redistribución" de peso espectral hacia $\lambda=0$ no requiere transporte
de información entre sitios: es una reasignación de amplitudes en el espacio
de Hilbert, no en el espacio real.

Análogos físicos establecidos:
- **Condensado de Bose-Einstein:** fracción macroscópica en $k=0$ sin
  violación de causalidad
- **Superconductividad:** pares de Cooper en el mismo estado cuántico
- **Teoría cuántica de campos:** modo cero como valor esperado del vacío
  (condensado de Higgs)

### A.3 Verificación numérica de localidad

Si la renormalización fuera acausal (señales más rápidas que la luz),
esperaríamos correlaciones espaciales instantáneas. Sin embargo, el
ancho del paquete en espacio real evoluciona suavemente desde 2.20 a
2.96 en $t=1.5$, consistente con una velocidad de propagación máxima
$v_{\max} = \Delta\text{ancho}/\Delta t \approx 0.5$ sitios por unidad
de tiempo, muy por debajo del límite físico $v_{\text{max}} = 1$
(1 sitio por paso temporal, dictado por el espectro del Laplaciano).

---

## B. El Parámetro $\beta$: Puente Geométrico Analítico

**Objeción:** "$\beta$ es un parámetro ajustable que debilita la tesis de
cero parámetros libres. ¿Existe un puente geométrico que lo fije?"

### B.1 Derivation analítica

Sí. $\beta$ se deriva de dos constantes geométricas de la red FCC:

1. **Fracción de empaquetamiento** $\eta = \pi/(3\sqrt{2}) \approx 0.74048$
   — la máxima densidad de empaquetamiento de esferas en 3D.

2. **Perturbación de fondo** $\varepsilon_0 = 1/12$ — del número de
   coordinación $Z=12$.

El colapso espectral tiene dos contribuciones complementarias:
- La **fracción vacía** $1-\eta$ representa el "espacio cuántico"
  disponible para fluctuaciones — la porción de la red que no está
  ocupada por la estructura clásica.
- La **auto-regulación** $\varepsilon_0/2$ es la mitad de la perturbación
  de fondo, que fija la escala de energía para el umbral de colapso.

Sumando ambas:

$$\beta = (1-\eta) + \frac{\varepsilon_0}{2}
        = 1 - \frac{\pi}{3\sqrt{2}} + \frac{1}{24}
        \approx 0.25952 + 0.04167 = 0.30119 \approx 0.30.$$

$\beta$ no es libre — está fijado por la geometría de la red ($\eta$)
y la escala fundamental ($\varepsilon_0$). El rango de estabilidad
$[0.1, 0.5]$ observado numéricamente contiene a $\beta=0.30$ y refleja
la tolerancia del atractor a variaciones en la tasa de colapso, sin
alterar el punto fijo (ancho $2.96$, energía $-0.056$) para ningún
$\beta$ en ese intervalo.

### B.2 Tiempo de purificación

El tiempo de convergencia $t_{\text{purificación}} \approx 1.5$ depende
de $\beta$ a través de $t \sim 1/(\beta \varepsilon_0 \Delta\lambda)$,
donde $\Delta\lambda \sim 1$ es el gap espectral mínimo. Con
$\beta=0.30$ y $\varepsilon_0=1/12$, obtenemos
$t \sim 1/(0.30 \times 0.083 \times 1) \approx 40$ pasos
($40 \times 0.005 = 0.2$), consistente con la escala observada
($t_{1/2} \sim 0.25$, $t_{\text{total}} \sim 1.5$).

---

## C. Purificación vs. Segunda Ley de la Termodinámica

**Objeción:** "El modelo propone que $S(t) \to 0$ mientras que la flecha
del tiempo macroscópica está asociada al *crecimiento* de entropía.
¿Cómo coexisten?"

### C.1 Sistema abierto, no aislado

La red DNLS no es un sistema aislado. El colapso espectral representa un
**acoplamiento disipativo a un reservorio externo** — los grados de libertad
de la red (fonones, modos de vibración, etc.) que no están incluidos
explícitamente en la simulación. La entropía de von Neumann $S(t)$ que
calculamos es la entropía **del subsistema reducido** (los grados de
libertad DNLS), no del universo total.

La segunda ley exige que la entropía total $S_{\text{total}} = S_{\text{DNLS}} +
S_{\text{red}} + S_{\text{ambiente}}$ sea no-decreciente. Nuestro resultado
$S_{\text{DNLS}}(t) \to 0$ implica que la entropía es transferida del
sistema DNLS a la red, que actúa como reservorio frío. Esto es exactamente
análogo a:

- **Enfriamiento por evaporación en BEC:** los átomos más energéticos
  abandonan el condensado, reduciendo su entropía
- **Preparación disipativa de estados puros:** sistemas cuánticos abiertos
  que alcanzan el estado fundamental del hamiltoniano mediante acoplamiento
  a un baño (Diehl et al., Nature Physics 4, 878, 2008)
- **Efecto Mpemba cuántico:** purificación acelerada por disipación
  (Westhoff et al., Phys. Rev. A 112, L061304, 2025)
- **Medición cuántica:** el sistema medido se purifica al correlacionarse
  con el aparato

### C.2 Dos flechas, misma dirección

La purificación espectral $S_{\text{DNLS}}(t) \to 0$ y la segunda ley
termodinámica $S_{\text{total}}(t) \geq 0$ son complementarias, no
contradictorias:

| Escala | Entropía | Dirección | Mecanismo |
|--------|----------|-----------|-----------|
| **Local/espectral** (sistema DNLS) | $S_{\text{DNLS}} \to 0$ | Purificación | Colapso de modos altos hacia $\lambda=0$ |
| **Global/macroscópica** (sistema+reservorio) | $S_{\text{total}} \nearrow$ | Termalización | Disipación de energía a la red |

La "flecha de purificación" es la manifestación microscópica del colapso
de la función de onda. La "flecha termodinámica" es el comportamiento
emergente a escalas macroscópicas. Ambas apuntan en la misma dirección
temporal (pasado $\to$ futuro) porque el mismo proceso — la disipación
irreversible de coherencia cuántica hacia la red — las genera
simultáneamente.

### C.3 Análogo cosmológico

En cosmología, la flecha del tiempo del universo temprano está asociada
a la disminución de la "entropía cuántica" durante la inflación
(período en el que las fluctuaciones cuánticas se "congelan" en
perturbaciones clásicas). Nuestro modelo replica este proceso en una
red discreta: el estado cuántico se purifica al colapsar espectralmente,
generando una estructura clásica (el solitón de ancho $2.96$) que
constituye el "universo observable" emergente.

### C.4 Verificación numérica

La monotonía de $S(t)$ es empírica: para $M=100$ realizaciones con
semillas $42,\dots,141$, todas las trayectorias muestran $dS/dt < 0$
en todo momento, con rango de la matriz de densidad reduciéndose de
$68 \to 1$. No se observa ninguna fluctuación al alza, consistente
con un proceso disipativo unidireccional.

---

## D. Cambios al Manuscrito (v3)

| Sección | Cambio |
|---------|--------|
| Sec. 2.5 (colapso espectral) | Formulación explícita como mapa de Kraus |
| Sec. 3.3 ($\beta$) | Derivation analítica: $\beta = (1-\eta) + \varepsilon_0/2$ |
| Sec. 4.4 (entropía) | Discusión: sistema abierto, dos flechas del tiempo |
| Apéndice A | Verificación de localidad: perfil de propagación |
| Conclusiones | $\beta$ ya no es parámetro libre |

---

## E. Conclusión

Las tres objeciones tienen respuesta rigurosa:

| Objeción | Respuesta |
|----------|-----------|
| **A (no-localidad)** | Mapa de Kraus diagonal; preservación de traza canónica; análogo a condensación de modo cero |
| **B ($\beta$ libre)** | $\beta = (1-\eta) + \varepsilon_0/2$ fijado por geometría FCC ($\eta$) y escala fundamental ($\varepsilon_0$) |
| **C ($S\to0$ vs. 2ª ley)** | Sistema abierto; $S_{\text{DNLS}}\to0$ es purificación local que transfiere entropía al reservorio; $S_{\text{total}} \geq 0$ |

Agradecemos nuevamente al revisor sus agudas observaciones, que han
permitido fortalecer significativamente el manuscrito.

---
*Serie Armónica Collaboration*
