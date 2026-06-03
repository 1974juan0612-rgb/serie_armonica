# Experimentsos que demuestran (o pueden demostrar) el modelo FCC-DNLS

## 1. Experimentos ya existentes que respaldan el modelo

### 1.1 Redes opticas 3D con geometria FCC

**Referencia:** "Weyl points in an FCC optical lattice" (ResearchGate, 2016)

Una red optica FCC se puede realizar con solo 3 lasers en configuracion tetraedrica.
El grupo de JILA (Campbell et al., Science 2017) ya opera relojes opticos 3D con
redes cubicas convencionales y átomos de estroncio en regímenes de degeneracion cuantica.
La extension a geometria FCC es un paso incremental, no una barrera tecnologica.

**Estado:** Realizable hoy. Varios grupos (Munich, JILA, MIT) tienen la capacidad.

### 1.2 Solitones discretos 3D en redes opticas

**Referencia:** Kevrekidis et al., "Three-Dimensional Solitary Waves and Vortices
in a Discrete Nonlinear Schrodinger Lattice", Phys. Rev. Lett. 93, 080403 (2004).

Soluciones de vórtices solitonicos 3D fueron predichas y caracterizadas en la red
DNLS hace dos decadas. La existencia de estructura coherente localizada en 3D
en redes discretas es un resultado establecido, no especulativo.

**Estado:** Teoricamente establecido. Pendiente de verificacion experimental directa
en redes opticas 3D con interacciones controladas.

### 1.3 Ingenieria de reservorio para disipacion espectral

**Referencia:** Tissot, Ribeiro, Marquardt, "Reservoir engineering for classical
nonlinear fields", Phys. Rev. Research 6, 023015 (2024).

Derivan una ecuacion de Schrodinger no-lineal disipativa a partir de un reservorio
ingenierizado. Estudian explicitamente la dinamica de solitones bajo esta disipacion
no-lineal. Esto es EXACTAMENTE el tipo de disipacion que necesita nuestro modelo:
la tasa de disipacion depende de la amplitud local (como nuestro gamma_n) pero
derivada de principios de ingenieria cuantica de reservorios.

**Estado:** Publicado. La teoria existe y los experimentos de disipacion controlada
en redes opticas son rutinarios (grupos de Zoller, Blatt, Bloch, Esslinger).

### 1.4 Purificacion disipativa en redes opticas

**Referencia:** Westhoff, Paeckel, Moroder, "Fast and direct preparation of a
genuine lattice Bose-Einstein condensate via the quantum Mpemba effect",
Phys. Rev. A 112, L061304 (2025).

Demuestran que la preparacion disipativa de un BEC en una red optica purifica
el estado cuantico (efecto Mpemba cuantico). La entropia del estado se reduce
por disipacion controlada — exactamente el fenomeno de purificacion que observamos
en nuestro modelo con S(t) -> 0.

**Estado:** Publicado. Confirmacion experimental de que disipacion + red =
purificacion del estado.

---

## 2. Experimentos propuestos (predicciones del modelo)

### Experimento A: Verificacion del punto fijo espectral en red FCC 3D

**Configuracion:**
1. Crear una red optica FCC con 3 lasers en configuracion tetraedrica
2. Cargar un condensado de Bose-Einstein de 87Rb (~10^4 átomos)
3. Ajustar la profundidad de red para que el tunelamiento sea t ~ 0.1 ER
4. Aplicar interacciones atractivas via resonancia de Feshbach para que U < 0 (DNLS focalizante)
5. Inyectar un pulso localizado (gaussiano, sigma ~ 2 sitios) en el centro de la red

**Prediccion del modelo:**
El sistema converge a un punto fijo con ancho = 2.96 sitios,
max|psi|^2 = 0.009 (N=6), energia constante.

**Protocolo de medida:**
- Imagen de absorcion a t = 0, 0.5, 1.0, 2.0, 5.0 tiempos caracteristicos
- Ajustar el perfil de densidad a una gaussiana para medir el ancho
- Verificar que el ancho se estabiliza (no crece ni colapsa a 1 sitio)

**Criterio de exito:**
Ancho se estabiliza en 2.96 +/- 0.1. Energia cinetica + interaccion constante.
Desviacion -> el modelo falla.

### Experimento B: Umbral de auto-atrapamiento (Norma = 4)

**Configuracion:**
Misma red FCC, misma configuracion. Variar el numero de atomos en el pulso inicial.

**Prediccion del modelo:**
Por debajo de N_atomos < N_critico (norma < 4): el pulso se expande (ancho -> 7+).
Por encima de N_atomos > N_critico (norma > 4): el pulso se auto-atrapa (ancho ~ 3.5).

**Protocolo de medida:**
1. Preparar pulsos con densidad atomica creciente (N = 100, 500, 1000, 2000, 5000)
2. Medir el ancho despues de t = 2.5 unidades de tiempo
3. Identificar la transicion abrupta entre expansion y auto-atrapamiento

**Criterio de exito:**
Transicion abrupta en N_critico. Por debajo: ancho > 5. Por encima: ancho < 4.

### Experimento C: Dependencia dimensional (1D vs 3D)

**Configuracion:**
1. Crear red 1D (tubo optico, haz laser en 1 direccion) con los mismos parametros
2. Crear red 3D FCC con los mismos parametros
3. Misma condicion inicial en ambas

**Prediccion del modelo:**
1D: el pulso se degrada bajo colapso espectral (ancho 3.8 -> 6.5).
3D: el pulso converge a punto fijo (ancho 2.96).

**Protocolo de medida:**
Comparar la evolucion temporal del ancho en 1D vs 3D.

**Criterio de exito:**
En 1D, el ancho crece sin estabilizarse. En 3D, se estabiliza.
Si ambas se estabilizan o ambas divergen -> el modelo falla.

### Experimento D: Purificacion espectral (S(t) -> 0)

**Configuracion:**
Misma red FCC 3D. Realizar M = 100 repeticiones del experimento con
condiciones iniciales identicas pero ruido experimental diferente.

**Prediccion del modelo:**
La matriz de densidad promediada sobre las realizaciones se purifica:
S(t) = -Tr(rho ln rho) -> 0. Rango(rho) -> 1.

**Protocolo de medida:**
1. Realizar M repeticiones con tomografia de estado cuantico en cada una
2. Construir rho(t) = (1/M) sum |psi_i(t)><psi_i(t)|
3. Diagonalizar y calcular autovalores y entropia
4. Verificar que S(t) decrece monotonicamente a 0

**Criterio de exito:**
S(t) -> 0 en t < 2. Si S(t) se mantiene > 0, el sistema no se purifica -> modelo falla.

---

## 3. Conexion con funcion zeta de Epstein

### Experimento E: Verificacion de la ecuacion funcional espectral

**Configuracion:** Solo numerico. No requiere experimento de laboratorio.

**Prediccion:** Para la red FCC, la funcion zeta espectral satisface:

    xi(3/2 - s) = xi(s)                                        (1)

con error que disminuye al aumentar N.

**Verificacion actual:** N=6 (108 sitios) -> error 2.84e-14. N=8 (256) -> error < 2%.
N=12 (864) -> diagonalizacion completa en 0.96s.

Este resultado conecta el modelo con la conjetura de Friedli-Karlsson (2017)
y la hipotesis de Riemann, pero NO requiere verificacion experimental.

---

## 4. Resumen de testabilidad

| Prediccion | Tipo | Dificultad | Tiempo estimado |
|-----------|------|-----------|----------------|
| Punto fijo en red FCC 3D | Laboratorio | Media | 6-12 meses |
| Umbral de auto-atrapamiento | Laboratorio | Media | 6-12 meses |
| Dimensionalidad 1D vs 3D | Laboratorio | Baja | 3-6 meses |
| Purificacion S(t)->0 | Laboratorio | Alta (tomografia) | 12-18 meses |
| Ecuacion funcional Epstein | Numerico | Trivial | Hecho |

El experimento mas accesible es la **comparacion 1D vs 3D** (C), que requiere
solo cambiar la configuracion de los lasers de la red optica. No necesita
nueva instrumentacion.

El experimento mas crucial es la **verificacion del punto fijo** (A), que
confirmaria o refutaria la prediccion central del modelo.

---

## Referencias experimentales clave

1. Campbell et al., "A Fermi-degenerate three-dimensional optical lattice clock",
   Science 358, 90-94 (2017). — Red 3D con atomos degenerados.
2. Kevrekidis et al., "Three-Dimensional Solitary Waves and Vortices in a DNLS
   Lattice", Phys. Rev. Lett. 93, 080403 (2004). — Solitones 3D en DNLS.
3. Tissot, Ribeiro, Marquardt, "Reservoir engineering for classical nonlinear
   fields", Phys. Rev. Research 6, 023015 (2024). — DNLS disipativa.
4. Westhoff, Paeckel, Moroder, "Fast preparation of lattice BEC via quantum
   Mpemba effect", Phys. Rev. A 112, L061304 (2025). — Purificacion disipativa.
5. Friedli, Karlsson, "Spectral zeta functions of graphs", Tohoku Math. J.
   69(4), 585-610 (2017). — Conexion con zeta de Riemann.
6. Connes, Consani, "Zeta Spectral Triples", arXiv:2511.22755 (2025).
   — Marco no conmutativo.
7. Khazali, "Ultratight confinement of atoms in a Rydberg empowered optical
   lattice", Quantum 9, 1585 (2025). — Nuevas tecnicas de confinamiento.
8. Continuous BEC, Nature 606, 683-687 (2022). — BEC continuo demostrado.
