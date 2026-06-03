# A: Galaxias masivas en el universo temprano (JWST)

**Pregunta:** Son las galaxias masivas a z>10 detectadas por JWST
un efecto de proyeccion espectral del soliton cosmologico FCC-DNLS?

**Hipotesis:** El corrimiento al rojo cosmologico no es expansion
metrica (FLRW) sino gradiente espectral dentro del perfil del
soliton. Objetos a distinto r tienen distinto corrimiento espectral
aparente.

---

## NOTA SOBRE LA CONVENCION DEL SIGNO DE z

En este modelo, el redshift se define como:
```
z(r) = (mean_lambda(observador) - mean_lambda(r)) / mean_lambda(r)
```
donde mean_lambda pondera cada modo espectral por su amplitud local.
z > 0 significa que el objeto tiene menor frecuencia espectral que
el observador (consistente con la convencion astronomica).

En A1b se observa z > 0 para r > 0 (periferia), con valores que van
de ~1.3 (r=1.25) a ~37 (r=3.25). Estos valores altos indican que la
frecuencia espectral cae rapidamente con la distancia al centro del
soliton, produciendo un gradiente espectral pronunciado.

**Consecuencia fisica:** Un objeto a r=3.25 emite con frecuencia
espectral ~37× menor que el observador. En FLRW esto corresponderia
a z=37, pero en nuestro modelo no hay expansion — es solo gradiente
de estado interno.

---

## Plan de trabajo

### Fase 1: Mapeo soliton -> metrica FLRW (teoria)

Derivar la metrica inducida por el perfil |psi(r)|^2 del soliton
FCC-DNLS y compararla con la metrica FLRW.

- [x] A1: Perfil radial del soliton y mapeo r->z
- [x] A1b: Redshift espectral dinamico durante el colapso
- [ ] Derivar g_{munu}^{soliton} como funcion de |psi(r)|^2 y xi
- [ ] Predecir desviaciones de H(z) respecto a LambdaCDM para z>6

### Fase 2: Simulacion de catalogos sinteticos (numerico)

Generar catalogos de galaxias sinteticos en el espacio de fase FCC
y comparar con datos JWST, Euclid, Roman.

- [x] A2: Replicacion morfologica espectral (z~2 -> z>10)
- [ ] Implementar distribucion de galaxias 3D en red FCC
- [ ] Aplicar proyeccion espectral (redshift como funcion de r)
- [ ] Generar diagramas H-R, funcion de masa, funcion de correlacion
- [ ] Comparar con datos JWST (PRIMER, CEERS, JADES)

### Fase 2.5: Confrontacion con datos observacionales

Validar las predicciones del modelo contra datos reales.

- [ ] Comparar H(z) predicha con cosmocronologia (relojes de galaxias)
- [ ] Contrastar con datos BAO (BOSS, DESI) para descartar modelo
- [ ] Verificar consistencia con abundancias de elementos ligeros (BBN)
- [ ] Calcular la edad del universo a z=10 en el modelo vs JWST

### Fase 3: Prediccion del eco espectral (observable)

Si el modelo es correcto, debe haber replicas periodicas en z de
estructuras cercanas, correspondientes a la periodicidad de la
supercelda N.

- [ ] Calcular el periodo espectral Delta_z ~ lambda_min / lambda_max
- [ ] Buscar en datos JWST galaxias con morfologias replicadas a
      intervalos regulares en z
- [ ] Proponer observacion de seguimiento con JWST Cycle 3

---

## Experimentos

| Archivo | Contenido | Estado |
|---------|-----------|--------|
| A1_perfil_soliton.py | Perfil radial y desviacion H(z) vs FLRW | Implementado |
| A1b_redshift_dinamico.py | Redshift espectral durante colapso (z>0 ahora) | Implementado |
| A2_replicacion_morfologica.py | Transformacion espectral de morfologias | Implementado |

## Hallazgos clave

1. **H(z) extremadamente baja:** desviacion de -75% a z=0.5, -93% a z=10
   respecto a FLRM. Implica edad del universo ~14× mayor a z=10 en el
   modelo. Consecuencia: resuelve el problema de galaxias masivas JWST,
   pero requiere verificacion con BAO y BBN.

2. **Redshift dinamico:** z(r,t) emerge del transitorio de colapso
   (t<0.5) y se congela en t~1.0. El universo observable estaria en
   un estado transitorio donde el colapso espectral aun no ha
   purificado todas las inhomogeneidades.

3. **Replicacion morfologica:** Galaxias a z>10 son transformaciones
   espectrales de galaxias a z~2. La similitud replicada es ~1.0
   para todas las morfologias (el filtro espectral es deterministico).

## Referencias

- Labbe et al. (2023) "A population of red candidate galaxies
  at z~12-20" Nature
- Finkelstein et al. (2023) "CEERS: Galaxies at z~8-15"
- Naidu et al. (2022) "JWST ERO: Massive galaxies at z~10"
- Boylan-Kolchin (2023) "Stress testing LambdaCDM with JWST"
