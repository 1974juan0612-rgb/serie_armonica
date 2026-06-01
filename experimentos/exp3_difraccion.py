import numpy as np
import requests

# =====================================================================
# EXP 3: MAPEO DE DIFRACCIÓN Y CANALES ELÁSTICOS (OCTATOPO)
# =====================================================================

print("=== EXPERIMENTO 3: BIFURCACIÓN DEL OCTATOPO Y DENSIDAD DE FLUJO ===")
print("Objetivo: Canalizar el tren de ondas armónicas por las 8 direcciones")
print("          de la junta elástica y medir la permeabilidad de la red FCC")
print("Restricción: Base 8 con acoplamiento FCC bajo carga Q²")
print("-" * 75)

# URL oficial de HepData (CERN) - Modo contingencia automático
url_hepdata = "https://www.hepdata.net/record/ins1358249?format=json"

print("-> Inyectando escala de presiones desde el CERN...")
try:
    response = requests.get(url_hepdata, timeout=10)
    data = response.json()
    valores_q2 = data['data_tables'][0]['independent_variables'][0]['values']
    puntos_q2 = [float(v['value']) for v in valores_q2[:10]]
    print("-> Presiones reales cargadas con éxito.")
except:
    print("-> [Aviso] Usando rampa de presión estandarizada de HERA (Fallback)...")
    puntos_q2 = [3.5, 5.0, 6.5, 8.5, 10.0, 12.0, 15.0, 20.0, 24.0, 32.0]

# =====================================================================
# 1. CONSTANTES DE LA ARQUITECTURA FCC — HEREDADAS DEL EXPERIMENTO 2
# =====================================================================
MÉTRICA_ÚNICA = 1.0
DIMENSION_LIMITE = 8.0                        # Cero trivial de Riemann
N_1 = 12                                      # Coordinación FCC
CAPAS_BUNKER = 30                             # Profundidad del búnker
tiempo_pi = np.pi
plano_100 = 100.0                             # Inercia superficial del plano (100)

# Empaquetamiento máximo de canicas (límite geométrico de compactación)
eta_max = np.pi / (3.0 * np.sqrt(2))          # ~0.740480

# Fricción base de la Caja 2
rho_base = abs(((N_1 / (N_1 + 1)) ** 2) - (1.0 / np.sqrt(2)))
Q_quarks = 3.0 / np.sqrt(N_1)

# =====================================================================
# 2. OCTATOPO — LAS 8 DIRECCIONES GEOMÉTRICAS DE LA JUNTA ELÁSTICA
# =====================================================================
# Vectores de los 8 canales de difracción en el plano (100) de la FCC
# Cada dirección corresponde a un vector de onda armónica en la base
octatopo = {
    "Eje X+":    np.array([ 1,  0,  0]),
    "Eje X-":    np.array([-1,  0,  0]),
    "Eje Y+":    np.array([ 0,  1,  0]),
    "Eje Y-":    np.array([ 0, -1,  0]),
    "Eje Z+":    np.array([ 0,  0,  1]),
    "Eje Z-":    np.array([ 0,  0, -1]),
    "Diagonal":  np.array([ 1,  1,  0]) / np.sqrt(2),
    "Antidiag":  np.array([-1,  1,  0]) / np.sqrt(2),
}

print(f"-> Octatopo desplegado: {len(octatopo)} canales de difracción activados")
print(f"-> Densidad de empaquetamiento FCC: {eta_max:.6f}")
print(f"-> Fricción angular de la red:      {rho_base:.6f}")
print(f"-> Inercia del plano (100):         {plano_100:.1f}")
print("-" * 75)

# =====================================================================
# 3. BUCLE PRINCIPAL — BARRIDO DE PRESIONES
# =====================================================================
print(f"{'Paso':<6}{'Q2 (GeV2)':<14}{'Permeab P':<14}{'Flujo B':<16}{'Fase Phi':<14}{'Armonica H':<14}{'Ocupac W':<14}")
print("-" * 75)

for idx, Q2 in enumerate(puntos_q2):
    # 3a. BIFURCACIÓN EN EL OCTATOPO
    # Cada dirección recibe un tren de ondas con fase modulada por Q²
    tren_ondas = {}
    for nombre_dir, vec_dir in octatopo.items():
        # Número de onda en esta dirección, modulado por la presión Q²
        k_dir = np.linalg.norm(vec_dir) * np.sqrt(Q2 + MÉTRICA_ÚNICA)
        # Fase acumulada en el recorrido por las 30 capas
        fase = 0.0
        for j in range(1, CAPAS_BUNKER + 1):
            atenuacion = (N_1 / (N_1 + j)) ** 2
            rho_capa = rho_base * (1.0 - (1.0 / (j + Q_quarks)))
            transmision = atenuacion * (1.0 - rho_capa)
            fase += transmision * np.sin(k_dir * j / CAPAS_BUNKER)
        tren_ondas[nombre_dir] = fase

    # 3b. MEDICIÓN DE PERMEABILIDAD (P)
    # Capacidad de la red FCC para transmitir vibración sincrónica
    # a través de las 30 capas, integrando la resistencia por fricción
    suma_transmision = 0.0
    for j in range(1, CAPAS_BUNKER + 1):
        nodos_capa = 10 * (j ** 2) + 2
        atenuacion_j = (N_1 / (N_1 + j)) ** 2
        rho_capa = rho_base * (1.0 - (1.0 / (j + Q_quarks)))
        transmision_capa = nodos_capa * atenuacion_j * (1.0 - rho_capa)
        # Factor de sincronía: qué tan alineada está la vibración con la red
        sincronia = np.cos(1.0 / (j + Q_quarks))
        suma_transmision += transmision_capa * sincronia

    # Permeabilidad: flujo que logra atravesar toda la red
    semi_eje = (DIMENSION_LIMITE / 2.0) / np.log10(10.0 + Q2)
    permeabilidad_P = (suma_transmision / CAPAS_BUNKER) / (eta_max * semi_eje)

    # 3c. TRAZADO DE DENSIDAD DE FLUJO (B)
    # Firma geométrica emergente regulada por la inercia del plano (100)
    # Se verifica si se forman pasillos de fuerza ordenados
    suma_flujo = 0.0
    for nombre_dir, fase_dir in tren_ondas.items():
        # Amplitud de la onda en esta dirección
        amplitud = abs(fase_dir)
        # Proyección sobre el plano de inercia (100)
        factor_plano = 1.0 + (amplitud / plano_100)
        suma_flujo += amplitud * factor_plano

    # Densidad de flujo emergente B
    flujo_B = (suma_flujo * eta_max) / (plano_100 * np.sqrt(Q2 + 1.0))

    # 3d. FASE GLOBAL — Coherencia del tren armónico
    fase_global = np.arctan2(
        sum(np.sin(f) for f in tren_ondas.values()),
        sum(np.cos(f) for f in tren_ondas.values())
    )

    # 3e. ARMÓNICA H — Número de onda del canal dominante
    canal_dominante = max(tren_ondas, key=lambda k: abs(tren_ondas[k]))
    armónica_H = abs(tren_ondas[canal_dominante]) * Q_quarks

    # 3f. OCUPACIÓN Ω — Fracción de canales con transmisión no nula
    canales_activos = sum(1 for f in tren_ondas.values() if abs(f) > 0.01)
    ocupacion_Omega = canales_activos / len(octatopo)

    print(f"{idx+1:<6}{Q2:<14.2f}{permeabilidad_P:<14.4f}{flujo_B:<16.4f}{fase_global:<14.4f}{armónica_H:<14.4f}{ocupacion_Omega:<14.2f}")

print("-" * 75)

# =====================================================================
# 4. LECTURAS PURAS DEL SILICIO — DIAGNÓSTICO FINAL
# =====================================================================

print("\n=== LECTURAS PURAS DEL SILICIO ===")
print("-" * 75)

# Mapa de coherencia direccional
coherencia_media = np.mean([abs(tren_ondas[d]) for d in octatopo])
print(f"Coherencia direccional media:           {coherencia_media:.6f}")

# Pasillos de fuerza ordenados (canales con flujo > umbral)
pasillos = [d for d in octatopo if abs(tren_ondas[d]) > coherencia_media]
print(f"Pasillos de fuerza ordenados:           {len(pasillos)}/{len(octatopo)}")
for p in pasillos:
    print(f"  -> Canal activo: {p:<12} | Amplitud: {abs(tren_ondas[p]):.4f}")

# Índice de jamming por presión
indice_jamming = 1.0 - (permeabilidad_P / (1.0 + permeabilidad_P))
print(f"Indice de jamming por presion:          {indice_jamming:.6f}")

# Verificación de simetría FCC (los 4 ejes principales deben ser equivalentes)
simetria_xy = abs(tren_ondas["Eje X+"]) - abs(tren_ondas["Eje Y+"])
simetria_yz = abs(tren_ondas["Eje Y+"]) - abs(tren_ondas["Eje Z+"])
print(f"Simetria FCC (D|X|-|Y|):                {simetria_xy:.6f}")
print(f"Simetria FCC (D|Y|-|Z|):                {simetria_yz:.6f}")

# Frecuencia natural del búnker
frecuencia_bunker = (DIMENSION_LIMITE * eta_max) / (N_1 * rho_base)
print(f"Frecuencia natural del búnker:          {frecuencia_bunker:.6f}")

print("-" * 75)
print("=== EXPERIMENTO 3 COMPLETO: DIFRACCIÓN Y CANALES ELÁSTICOS ===")
