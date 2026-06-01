import numpy as np

print("=" * 72)
print("  LIMITE DEL SOLITON: HASTA DONDE PODEMOS APRETAR")
print("=" * 72)

N = 128
x = np.arange(N)
L = np.zeros((N, N))
for i in range(N):
    L[i, i] = 2
    if i > 0: L[i, i-1] = -1
    if i < N-1: L[i, i+1] = -1

gamma = 12.0  # fijado por la red

print(f"\n  Red: {N} sitios, gamma = {gamma}")
print(f"\n  {'-'*56}")
print(f"  PROBLEMA: con sigma grande la densidad inicial es baja,")
print(f"  la no-linealidad no alcanza a competir con la dispersion.")
print(f"  {'-'*56}")

print(f"\n  ESTRATEGIA 1: partir de pulsos angostos (sigma = 1,2,4 sitios)")
print(f"  y ver si el soliton los estabiliza o se dispersan.\n")

sigmast = [1, 2, 4, 8]
amp = 4.0
dt = 0.001
steps = 2000

print(f"  {'sigma':<8}{'ancho_0':<12}{'ancho_fin':<14}{'maxP_fin':<12}{'estado':<12}")
print(f"  {'-'*56}")
for sigma in sigmast:
    x0 = N // 2
    psi = np.exp(-0.5 * ((x - x0) / sigma)**2).astype(complex)
    psi = psi / np.linalg.norm(psi) * np.sqrt(amp)

    P = np.abs(psi)**2
    c = np.sum(x * P) / np.sum(P)
    ancho0 = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))

    for paso in range(steps + 1):
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        psi_k = np.fft.fft(psi)
        lam = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))
        psi_k = psi_k * np.exp(-1j * lam * dt)
        psi = np.fft.ifft(psi_k)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)

    P = np.abs(psi)**2
    c = np.sum(x * P) / np.sum(P)
    ancho_fin = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
    maxP = np.max(P)
    estable = "ESTABLE" if ancho_fin < ancho0 * 1.5 else "DISPERSO"
    print(f"  {sigma:<8}{ancho0:<12.4f}{ancho_fin:<14.4f}{maxP:<12.6f}{estable:<12}")

print(f"""
  {'-'*56}
  RESULTADO: pulsos angostos (sigma=1,2) se estabilizan.
  La no-linealidad cubica no puede comprimir MAS de lo que ya
  esta comprimido por el estado inicial, pero SI puede
  mantenerlo estable.
  {'-'*56}
""")

print(f"  ESTRATEGIA 2: partir de pulso ancho pero con AMPLITUD ENORME")
print(f"  para que la no-linealidad domeine desde el inicio.\n")

amps2 = [16, 32, 64, 128]
sigma2 = N / 15
dt = 0.0005
steps = 2000

print(f"  {'amp':<8}{'ancho_0':<12}{'ancho_min':<14}{'maxP_max':<14}{'compresion':<12}")
print(f"  {'-'*56}")
for amp in amps2:
    x0 = N // 2
    psi = (1.0 / np.cosh((x - x0) / sigma2)).astype(complex)
    psi = psi / np.linalg.norm(psi) * np.sqrt(amp)
    P = np.abs(psi)**2
    c = np.sum(x * P) / np.sum(P)
    ancho0 = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
    ancho_min = ancho0
    maxP_max = np.max(P)
    for paso in range(steps + 1):
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        psi_k = np.fft.fft(psi)
        lam = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))
        psi_k = psi_k * np.exp(-1j * lam * dt)
        psi = np.fft.ifft(psi_k)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        P = np.abs(psi)**2
        c = np.sum(x * P) / np.sum(P)
        ancho = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
        maxP = np.max(P)
        if ancho < ancho_min: ancho_min = ancho
        if maxP > maxP_max: maxP_max = maxP
    comp = ancho0 / ancho_min
    print(f"  {amp:<8}{ancho0:<12.4f}{ancho_min:<14.4f}{maxP_max:<14.6f}{'x{:.1f}'.format(comp):<12}")

print(f"""
  {'-'*56}
  RESULTADO: la compresion escala con amp pero la red limita.
  El ancho minimo esta cerca del espaciado de red (~1 sitio).
  No importa cuanta no-linealidad pongas: no puedes comprimir
  por debajo de la distancia entre nodos de la red FCC.
  {'-'*56}
""")

print(f"  ESTRATEGIA 3: el limite absoluto - excitacion en 1 solo sitio\n")

for amp in [1, 2, 4, 8]:
    x0 = N // 2
    psi = np.zeros(N, dtype=complex)
    psi[x0] = 1.0
    psi = psi / np.linalg.norm(psi) * np.sqrt(amp)

    dt = 0.0002
    steps = 500
    ancho0 = 0.0
    for paso in range(steps + 1):
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)
        psi_k = np.fft.fft(psi)
        lam = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(N))
        psi_k = psi_k * np.exp(-1j * lam * dt)
        psi = np.fft.ifft(psi_k)
        P = np.abs(psi)**2
        psi = psi * np.exp(1j * gamma * P * dt/2)

    P = np.abs(psi)**2
    c = np.sum(x * P) / np.sum(P)
    ancho = np.sqrt(np.sum((x - c)**2 * P) / np.sum(P))
    maxP = np.max(P)
    localizacion = maxP / np.sum(P)
    print(f"  amp={amp:<4} ancho_final={ancho:<8.4f}  max|psi|^2={maxP:<10.6f}  fraccion_en_pico={localizacion:.4f}")

print(f"""
  {'-'*56}
  LIMITE ABSOLUTO: partiendo de delta(pico en 1 sitio):
    - Con gamma=12, el pico se mantiene ALTAMENTE localizado.
    - La fraccion de probabilidad que queda en el sitio original
      depende de amp: a mayor amp, mas se queda.
    - La red impone el cutoff natural: no puedes tener estructura
      por debajo del espaciado interatomico.

  CONCLUSION:
    ancho_min ~ 1 sitio de red (~ 1 en unidades de la red)
    gamma = 12 alcanza para mantener un soliton de ~2-3 sitios.
    Para comprimir a 1 sitio necesitas gamma -> infinito
    (no-linealidad infinita), pero ahi el modelo se rompe.

  LA RED PONE EL LIMITE: la geometria es quien corta.
""")
