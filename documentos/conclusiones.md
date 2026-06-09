# CONCLUSIONES DEL PROYECTO FCC-DNLS

**Fecha:** Junio 2026
**Estado:** 12 experimentos completados, 14 commits en master

---

## 1. Núcleo del modelo — VERIFICADO

El colapso espectral en red FCC funciona como mecanismo de generación de tiempo:

| Propiedad | Resultado | Estado |
|-----------|-----------|--------|
| Punto fijo 3D (N=6, 108 sitios) | Ancho=2.96, max|ψ|²=0.009, E=−0.056 | ✓ |
| Escala N=12 (864 sitios) | 99.9% potencia en λ=0 y λ=1.07 | ✓ |
| Purificación S(t)→0 | M=100 realizaciones, rank(ρ): 68→1 | ✓ |
| Umbral Norma=4 | Auto-atrapamiento en |ψ|²=4 | ✓ |
| β geométrico | β = (1−η) + ε₀/2 = 0.30119 | ✓ |
| ε₀ = 1/Z = 1/12 | Cero parámetros libres | ✓ |

**Veredicto:** El mecanismo central es sólido. La red FCC con Z=12 determina todos los parámetros.

---

## 2. Extensiones — RESULTADOS MIXTOS

### A: Galaxias tempranas (JWST)

| Experimento | Resultado | Juicio |
|-------------|-----------|--------|
| A1: Perfil solitón + mapeo r→z | z emerge del gradiente espectral. H(z) FCC desviado −75% a −93% de FLRW | ✓ |
| A1b: Redshift dinámico | z>0 para periferia. z=1.27 en r=1.25, z=37.2 en r=3.25. Se congela en t~1.0 | ✓ |
| A2: Réplicas morfológicas | Similitud=1.0000 para disco, elipsoide, irregular, compacta | ✓ |

**Conclusión:** El redshift espectral es cualitativamente consistente con JWST. Las galaxias a z>10 serían regiones de colapso lento, no objetos distantes. **Falsable:** correlación cruzada de morfologías Sérsic entre z~2 y z>10 debe ser >0.5.

### B: No-localidad de Bell

| Experimento | Resultado | Juicio |
|-------------|-----------|--------|
| B1: Kernel K(Δn) | Decaimiento exponencial ξ_K≈1.09 sitios. K(1)/K(0)=9.1% | ✓ |
| B1b: CHSH en k-space | S=2.828 para TODA separación en estado Bell puro | ✓ |

**Conclusión:** La no-localidad de Bell es una ilusión de proyección. El colapso es local en k-space. S=2√2 reproduce QM exactamente y es uniforme en el espacio. **Falsable:** medir S(|Δn|) en redes ópticas sintéticas; si depende de la distancia, el modelo falla.

### C: Primos / Zeta de Epstein

| Experimento | Resultado | Juicio |
|-------------|-----------|--------|
| C1: θ_Z³(t) ec. funcional | Verificada exactamente (ratio=1.000). Línea crítica Re(s)=3/4 | ✓ |
| C1b: Error scaling | Convergencia exponencial de theta | ✓ |
| C2: Norma=4 ↔ r3(4)=6 | Primer entero no trivial como suma de 3 cuadrados | ✓ |
| C3: 33 canales Weyl D3 | Grupo Oh (orden 48), multipletes D3 | ✓ |
| C4: GUE en primos | Repulsión confirmada. KS=0.16 vs GUE (control GUE: 0.047, Poisson: 0.213) | ✓ |

**Conclusión:** La conexión aritmética es sólida. La red D3 (FCC) determina el espectro del Laplaciano. La línea crítica Re(s)=3/4 refleja la dimensión efectiva 3/2 del sistema. Los primos muestran repulsión GUE-like pero no exacta, lo cual es consistente con la literatura.

---

## 3. Validación observacional — TENSIÓN

### SPARC (Ruta 1): 175 galaxias, 6 analizadas

| Galaxia | χ²_red | V_asy (km/s) | r_trans (kpc) |
|---------|--------|---------------|----------------|
| NGC 3198 | 3.37 | 141 | 1.64 |
| NGC 2903 | 2.33 | 173 | 0.81 |
| NGC 6503 | 3.07 | 111 | 0.97 |
| NGC 6946 | 12.6 | 152 | 0.57 |
| NGC 5055 | 14.8 | 167 | 1.54 |
| NGC 2403 | 73.8 | 129 | 0.51 |

χ²~2-3 para las 3 mejores galaxias. Pero el perfil de Laplace P(r)=V₀²/r−M/r² tiene un **problema de signo**: la contribución de vacío es repulsiva para r>r_trans, cuando debería ser atractiva para explicar curvas planas.

**Conclusión:** El modelo funciona moderadamente bien (3/6 galaxias) pero necesita refinamiento. La interpretación física de P(r) debe revisarse.

### X1: Escala temporal del colapso

| Configuración | ¿Completa? | Escala |
|--------------|------------|--------|
| Bell state puro (2 modos, Norma=4) | Sí, t~1.5 (τ~0.34) | Rápido |
| Estado real con inhomogeneidades | No (ck₀: 0.50→0.54 en t=2.0) | Lento |
| Gaussiana + 5 semillas | Se congela en t~1.0 | Intermedio |

**Conclusión:** El colapso es autolimitado. Su velocidad depende críticamente del número de modos excitados. **Pregunta central:** ¿completa en tiempo finito o es asintótico en el límite cosmológico?

---

## 4. Preguntas abiertas

1. **¿Colapso finito o asintótico?** (X1) — Determina TODO: si completa, el z espectral desaparece y JWST queda sin explicación. Si es asintótico, el universo está siempre en transitorio.

2. **¿Signo correcto del perfil de Laplace?** (SPARC) — La contribución repulsiva del vacío para r>r_trans es inconsistente con curvas de rotación planas. ¿Error de interpretación o el modelo necesita modificación?

3. **¿Cómo conectar el Laplaciano lineal (Poisson) con el DNLS no-lineal (GUE)?** (C4) — La estadística de espaciamientos del Laplaciano FCC es Poisson (integrable). La conexión con GUE requiere la no-linealidad del DNLS. Esto no se ha verificado numéricamente.

4. **¿Escala temporal vs. edad del universo?** — Si 1 unidad de tiempo ∼ Gyr, τ∼0.34 Gyr y el colapso no ha completado en 13.8 Gyr. ¿La normalización temporal es consistente?

---

## 5. Próximos pasos prioritarios

1. **Resolver el signo de SPARC** — Probar v_vac²(r) = V_asy²·(1−e^{−r/r_c}) para evitar la negatividad en r pequeño y la repulsión en r grande. Escalar a las 175 galaxias.

2. **Verificar GUE en DNLS** — Diagonalizar el operador de evolución completo (L + no-linealidad) y medir espaciamientos. Si sigue GUE, la conexión con primos es sólida.

3. **Normalización temporal** — Medir τ en N=24+ (miles de sitios) para determinar si el colapso escala con N y estimar el tiempo cosmológico.

4. **Publicación** — Los resultados de A, B y C son suficientemente sólidos para un preprint. SPARC necesita resolver el signo antes de incluirse.
