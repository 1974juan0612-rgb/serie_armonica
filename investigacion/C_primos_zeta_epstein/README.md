# C: Numeros primos, funcion zeta de Epstein y red D3 (FCC)

**Pregunta:** Por que la ecuacion funcional de Epstein para la red
D3 (isomorfa a FCC) tiene linea critica en Re(s) = 3/4, y que
relacion tiene con el umbral de auto-atrapamiento Norma=4?

**Hipotesis:** La red D3 tiene la funcion zeta de Epstein:
  zeta_D3(s) = sum_{(n1,n2,n3) in Z^3 \ {0}} w(n) / (n1^2+n2^2+n3^2)^s
con pesos w(n) = 1 si n1+n2+n3 es par (FCC).

La ecuacion funcional es:
  xi(s) = pi^{-s} * Gamma(s) * zeta_D3(s)
  xi(3/2 - s) = xi(s)

Esto implica que los ceros no triviales estan sobre Re(s) = 3/4.
El umbral Norma=4 del FCC-DNLS corresponde al primer autovalor del
Laplaciano FCC con degeneracion aritmeticamente distinguida.

---

## Plan de trabajo

### Fase 1: Ecuacion funcional de Epstein

- [ ] C1: Verificar ecuacion funcional numericamente
- [ ] C1b: Escalamiento del error con N: O(N^{-1/2})
- [ ] Encontrar los primeros ceros de xi(s) en Re(s) = 3/4

### Fase 2: Umbral aritmetico

- [ ] C2: Relacion entre Norma=4 y representaciones de sumas de 3 cuadrados
- [ ] Mostrar que lambda=4 tiene degeneracion 6 (primera representacion no trivial)
- [ ] Conectar con el teorema de los tres cuadrados de Legendre-Gauss

### Fase 3: Multipletes de Weyl

- [ ] C3: Clasificar los primeros 33 modos por representaciones del grupo de Weyl D3
- [ ] Mostrar estructura de multipletes: 1, 3, 6, 8, 12, ...
- [ ] Relacion con coeficientes de Clebsch-Gordan

---

## Experimentos

| Archivo | Contenido | Estado |
|---------|-----------|--------|
| C1_ecuacion_funcional.py | Verificacion ecuacion funcional Epstein | Pendiente |
| C1b_error_scaling.py | Escala O(N^{-1/2}) del error | Pendiente |
| C2_umbral_aritmetico.py | Norma=4 y sumas de 3 cuadrados | Pendiente |
| C3_multipletes_Weyl.py | Multipletes D3 y 33 canales | Pendiente |

## Referencias

- Epstein (1903) "Zur Theorie allgemeiner Zetafunktionen"
- Friedli & Karlsson (2017) "Zeta functions for lattices"
- Terras (1985) "Harmonic Analysis on Symmetric Spaces"
- Conway & Sloane (1999) "Sphere Packings, Lattices and Groups"
- Weil (1967) "Basic Number Theory"
