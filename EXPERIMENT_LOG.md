# EXPERIMENTO SERIE ARMONICA

## Estado actual — donde lo dejamos

### Lo que funciona (validado)
- Red FCC: Laplaciano exacto, autovalores en [0, 16], 12 vecinos/sitio
- Espectro coincide con formula analitica (lambda(k) = 12 - 4[cos cos + cos cos + cos cos])
- epsilon_0 = 1/12: unico valor que balancea colapso del centro (~0) con limpieza de bordes (~1)
- gamma = 1/epsilon_0 = 12: unico valor que cierra el modelo (cero parametros libres)
- Soliton DNLS desde pulso localizado (ancho 3.8 -> 2.7 con gamma=12)
- Auto-arranque: gamma_n = eps0/(eps0+|psi|^2) funciona sin observador externo

### La barrera (no resuelto)
- Colapso aleatorio + DNLS = caos (se destruyen mutuamente)
- Colapso filtrante = sobre-enfoque (todo a 1 sitio, no es soliton)
- La integracion de colapso (tiempo) y no-linealidad (estructura) no se logro
- Esto reproduce el problema de la medida: lo cuantico y lo clasico no se soldan

### Archivos clave
experimentos/
  exp_fcc_fundamento.py        - FCC Laplaciano, autovalores, zeta espectral
  exp_autoarranque.py          - Self-starting con epsilon_0
  exp_soliton_red.py           - DNLS soliton (corregido, gamma=4.0)
  exp_unificado.py             - Modelo completo con explicacion matematica
  exp_validacion.py            - 4 validaciones (2 pasan, 2 fallan)
  exp_colapso_focalizante.py   - Intento de colapso constructivo
  exp_limite_soliton.py        - Limite de compresion del soliton
  exp_prueba_tiempo.py         - Pruebas de tiempo/inercia
  exp_tiempo_inercia.py        - Original: tiempo de inercia

### Proximo paso (pendiente)
- Separar modos: colapso solo en frecuencias altas (ruido), preservar bajas (soliton)
- Probar si asi se salva la barrera
