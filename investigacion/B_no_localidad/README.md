# B: No-localidad de Bell como ilusion de k-espacio

**Pregunta:** Es la no-localidad cuantica un artefacto de proyectar
sobre la base real el colapso espectral (que ocurre en k-espacio)?

**Hipotesis:** El colapso ocurre en un unico modo k_0 en el espacio
de Fourier. La correlacion instantanea entre detectores separados
aparece porque la transformada de Fourier del kernel K(Delta n) es
no-local en x-espacio. Si midieramos directamente en k-espacio
(espectroscopia de momentos), veriamos variables locales clasicas.

---

## NOTA CONCEPTUAL

En el modelo FCC-DNLS:
- El colapso es LOCAL en k-espacio: actua modo a modo con tasa
  Gamma_k = sigma(lambda_k) * eps0 / (eps0 + |c_k|^2).
- En x-espacio, el mismo proceso se describe por un kernel:
  K(r, r') = Sigma_k e_k(r) * Gamma_k * e_k*(r')
- K(Delta n) tiene SOPORTE COMPACTO si Gamma_k es suave en k.
- Dos detectores en r1, r2 miden en base real. La correlacion
  instantanea aparece porque ambos proyectan sobre el mismo modo k_0,
  y K(r1, r') y K(r2, r') tienen overlap no nulo.

**Prediccion falsable:** Si el parametro CHSH S(theta1, theta2) se
mide en funcion de la separacion |r1 - r2|, debe decaer desde
S ~ 2.828 (violacion maxima) para |Delta n| < xi (tamano soliton)
hasta S < 2 (no violacion) para |Delta n| >> xi.

---

## Plan de trabajo

### Fase 1: Kernel de colapso espectral

- [ ] B1: Calcular K(Delta n) para FCC 3D. Verificar soporte compacto.
- [ ] Caracterizar decaimiento: |K(Delta n)| ~ exp(-|Delta n|/xi_K)
- [ ] Comparar xi_K con xi_soliton (ancho del perfil)

### Fase 2: Violacion CHSH simulada

- [ ] B1b: Simular pares de fotones correlacionados en modo k_0
- [ ] Calcular S(theta1, theta2) en funcion de |r1 - r2|
- [ ] Mostrar transicion S > 2 -> S < 2 para |Delta n| > xi

### Fase 3: Conexion con experimentos

- [ ] Comparar con Tissot et al. (PRR 6, 023015, 2024) - redes FCC
- [ ] Proponer implementacion en redes opticas sinteticas

---

## Experimentos

| Archivo | Contenido | Estado |
|---------|-----------|--------|
| B1_kernel_colapso.py | Kernel K(Delta n) en FCC 3D | Implementado |
| B1b_violacion_CHSH.py | Parametro S vs separacion | Pendiente |

## Referencias

- Tissot et al. (2024) "Local description of entanglement through
  synthetic non-local interactions" PRR 6, 023015
- Bell (1964) "On the Einstein-Podolsky-Rosen paradox" Physics 1, 195
- Aspect et al. (1982) "Experimental test of Bell inequalities
  using time-varying analyzers" PRL 49, 1804
- Brunner et al. (2014) "Bell nonlocality" RMP 86, 419
