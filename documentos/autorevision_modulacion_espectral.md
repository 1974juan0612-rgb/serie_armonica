# AUTO-REVISION: Modulacion Espectral como Canal de Comunicacion

**Manuscrito:** FCC-DNLS v3 (colapso espectral, kernel K, modulacion)
**Autor:** Serie Armonica Collaboration (auto-revision)
**Fecha:** Junio 2026

---

## Resumen Ejecutivo

El experimento `exp_modulacion_espectral.py` demuestra que una
perturbacion local en x_1 altera globalmente los coeficientes
espectrales |c_k|^2, y que el kernel de colapso K en x_2 cambia
230% instantaneamente, independientemente de la distancia
(x_1-x_2 = 6.93 sitios). A primera vista, esto sugiere un canal
de comunicacion espectral.

**Esta auto-revision identifica 7 objeciones que debilitan o
matizan la afirmacion.** El fenomeno es real, pero su utilidad
como canal de comunicacion tiene limites estrictos que no fueron
explorados en el experimento original.

---

## Objecion 1: Violacion del Teorema de No-Comunicacion

**El problema.** En mecanica cuantica estandar, el teorema de
no-comunicacion (Eberhard 1978, Ghirardi 1980) establece que
operaciones locales en A no cambian la matriz de densidad reducida
en B. Nuestro modelo parece violarlo.

**Auto-critica.** El teorema de no-comunicacion asume:
1. Evolucion unitaria entre mediciones
2. Mediciones projectivas locales (POVMs)
3. Comunicacion clasica como unico canal lateral

Ninguna de estas condiciones se cumple en el modelo FCC-DNLS:
- La evolucion NO es unitaria (el colapso espectral es disipativo)
- No hay mediciones projectivas; la dinamica es continua y clasica
  (|psi(x,t)|^2 es un campo clasico legible)
- El colapso actua GLOBALMENTE en la base espectral

**Veredicto.** El modelo NO viola el teorema de no-comunicacion
porque opera en un regimen fisico diferente (sistema abierto
disipativo con dinamica clasica-cuanticas hbrida). La prediccion
de un canal espectral es NUEVA FISICA, no un error.
Pero esto mismo la hace sospechosa: predice algo que la fisica
estandar excluye. Necesita validacion experimental independiente.

---

## Objecion 2: El Factor de Perturbacion 2.0 es Irrealista

**El problema.** El experimento multiplica por 2 la amplitud en
x_1. Esto requiere inyectar energia:

```
Delta E = gamma/2 * (|psi_pert(x_1)|^4 - |psi_base(x_1)|^4)
        = 6 * (0.036^2 - 0.0093^2)
        = 6 * (0.00130 - 0.000086)
        = 0.0073 unidades
```

Una perturbacion realista (10% de modulacion) daria:

```
Delta E ~ 6 * ((1.1*0.0093)^4 - 0.0093^4) ~ 6 * 0.0000013 ~ 8e-6
```

y el cambio en |psi(x_2)|^2 seria ~ 2.5e-6, comparable al ruido
numerico (1e-16 de precision de maquina).

**Veredicto.** El canal solo funciona con modulaciones grandes
(factor > 1.5). Para modulaciones pequenas (< 10%), la senial se
pierde en el ruido de redondeo. La capacidad del canal es
proporcional a la energia que Alice puede inyectar.

---

## Objecion 3: El Ruido Estocastico del Colapso No Fue模型izado

**El problema.** En Fase 4, todas las realizaciones base dan
|psi(x_2)|^2 = 0.00925926 exactamente. Esto significa que el
punto fijo es DETERMINISTA: el ruido estocastico (theta aleatorio)
no afecta las amplitudes estacionarias.

**Auto-critica.** Esto es correcto para el punto fijo, pero el
experimento no modelo el escenario realista de comunicacion:
Alice modula su amplitud en t=0, y Bob mide en t=delta_t,
ANTES de que el sistema alcance el nuevo punto fijo.

Durante el transitorio (0 < t < 1.0), el ruido estocastico SI
afecta la trayectoria. El experimento solo midio en t=1.0
(200 pasos), cuando el transitorio ya termino. Para comunicacion
en tiempo real, Bob no puede esperar a que el sistema se
estabilice.

**Veredicto.** El experimento subestima el ruido al medir solo
en el punto fijo. La capacidad del canal en regimen transitorio
es menor que la reportada.

---

## Objecion 4: El Observable es un Campo Clasico, no Cuantico

**El problema.** En el modelo, |psi(x,t)|^2 es un campo clasico
que Bob puede leer continuamente sin colapsar el estado. Esto
NO es una medicion cuantica - es la lectura de un campo clasico.

**Implicacion.** El "canal de comunicacion espectral" es
equivalente a un canal clasico con acoplamiento global. La
senal viaja instantaneamente en el sentido de que el kernel K
se reconfigura globalmente, pero Bob solo puede leerla
localmente con una tasa limitada por el tiempo de integracion.

**Analogia exacta.** Es como tener un fluido incompresible en un
recipiente: si Alice empuja en un punto, el nivel cambia en todo
el recipiente instantaneamente (en la aproximacion de fluido
incompresible). Pero la senal que Bob detecta es el cambio de
nivel, que requiere tiempo para estabilizarse.

**Veredicto.** El canal es clasico, no cuantico. No hay ventaja
fundamental sobre un canal clasico tradicional, excepto la
geometria global de la conexion.

---

## Objecion 5: Dependencia de la Condicion de Auto-atrapamiento

**El problema.** El experimento usa norma = 4 (auto-atrapado).
En el regimen de expansion (norma < 4), el punto fijo no existe
y la dinamica es caotica. El canal espectral solo funciona en
el regimen de auto-atrapamiento.

**Veredicto.** El canal depende de la existencia del atractor
global. Si el sistema no esta auto-atrapado, la senial se
dispersa. Esto limita el canal a sistemas con amplitud
suficiente (norma >= 4).

---

## Objecion 6: Problema de Sincronizacion

**El problema.** El experimento asume que Alice y Bob comparten
un reloj comun (saben que t=0 es cuando Alice perturba). En un
escenario real:

1. Bob no sabe cuando Alice comienza a modular
2. Bob no sabe cuanto duro la modulacion (1 bit o 0 bits?)
3. La modulacion de Alice cambia el punto fijo, pero Bob no
   tiene una "linea base" de referencia

**Solucion teorica.** Alice y Bob necesitarian un protocolo de
sincronizacion, lo que requiere un canal clasico auxiliar.
Esto reintroduce la limitacion de velocidad de la luz para la
sincronizacion, anulando la ventaja de la modulacion espectral.

**Veredicto.** El canal espectral requiere sincronizacion
externa, lo que limita su aplicacion practica para comunicacion
superluminical. Sin embargo, para comunicacion dentro del
volumen del soliton (distancia < 3 sitios), la sincronizacion
podria ser inherente.

---

## Objecion 7: El Kernel K Depende del Estado, No es un Medio Fijo

**El problema.** El kernel K(Delta n) depende de |c_k|^2, que
cambia con cada modulacion. Esto significa que el canal NO es
lineal ni estacionario: la respuesta en x_2 depende del estado
global del sistema, que Alice esta modulando.

**Implicacion.** Si Alice modula para enviar un bit, cambia el
kernel. Si luego modula para enviar otro bit, el nuevo kernel
es diferente. La comunicacion requiere un modelo del canal que
evoluciona con el propio mensaje - una no-linealidad que puede
generar interferencia entre bits sucesivos (memoria de canal).

**Veredicto.** El canal es no-lineal y con memoria. La capacidad
de Shannon de un canal asi es dificil de calcular y puede ser
significativamente menor que la de un canal lineal sin memoria.

---

## Tabla Resumen

| Objecion | Severidad | Estado |
|----------|-----------|--------|
| 1. No-comunicacion | Alta | Refutada (regimen distinto) |
| 2. Factor 2 irrealista | Alta | Confirmada: SNR cae con modulacion pequena |
| 3. Ruido transitorio | Media | No modelado en el experimento |
| 4. Campo clasico | Baja | Es una caracteristica, no un bug |
| 5. Dependencia del umbral | Media | Canal solo funciona para norma >= 4 |
| 6. Sincronizacion | Alta | Requiere canal clasico auxiliar |
| 7. No-linealidad/memoria | Alta | Canal no estacionario |

---

## Conclusion de la Auto-Revision

**El fenomeno de modulacion espectral es REAL** - la perturbacion
en x_1 altera globalmente |c_k|^2 y el kernel K en x_2. El
cambio de 230% en gamma local en x_2 es innegable.

**Sin embargo, como canal de comunicacion, las limitaciones son
severas:**

1. Solo funciona en regimen de auto-atrapamiento (norma >= 4)
2. Requiere modulaciones grandes (factor > 1.5) para SNR > 1
3. El ruido transitorio limita el bit rate real
4. La sincronizacion requiere un canal clasico auxiliar
5. El canal es no-lineal y con memoria

**Utilidad cientifica real:**
- NO sirve para comunicacion superluminical practica
- SI sirve como laboratorio teorico para entender la relacion
  entre localidad espectral y no-localidad espacial
- SI predice un efecto medible en redes opticas FCC: la
  correlacion de largo alcance en el ruido de colapso
- SI abre la pregunta de si la no-localidad cuantica de Bell
  es un caso particular de este fenomeno espectral

**El "pajaro vuela", pero con alas de plomo.** El resultado
principal no es la comunicacion, sino la demostracion de que
la no-localidad aparente en x-espacio es localidad en k-espacio.
Eso es suficiente para una publicacion en Physical Review D.

---
*Serie Armonica Collaboration (auto-revision critica)*
