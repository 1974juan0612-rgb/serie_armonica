"""
SPARC ajuste completo 175 galaxias con perfil exponencial.
Modelo: v_obs^2 = Y*(v_disk^2 + v_bul^2) + v_gas^2 + V_asy^2 * (1 - exp(-r / r_c))
El perfil exponencial evita el problema de signo del perfil de Laplace.
"""
import numpy as np, math, sys, time
from scipy.optimize import curve_fit

data_raw = open(r"C:\Users\famil\Desktop\serie_armonica\SPARC\table2.dat", "rb").read().decode("latin-1")
lines = [l for l in data_raw.split("\n") if l.strip() and len(l) >= 76]

galaxies = {}
for line in lines:
    try:
        name = line[0:11].strip()
        dist = float(line[12:18].strip())
        rad  = float(line[19:25].strip())
        vobs = float(line[26:32].strip())
        evobs = float(line[33:38].strip())
        vgas = float(line[39:45].strip())
        vdisk = float(line[46:52].strip())
        vbul = float(line[53:59].strip())
    except ValueError:
        continue
    if name not in galaxies:
        galaxies[name] = {"dist": dist, "R": [], "Vobs": [], "eVobs": [], "Vgas": [], "Vdisk": [], "Vbul": []}
    galaxies[name]["R"].append(rad)
    galaxies[name]["Vobs"].append(vobs)
    galaxies[name]["eVobs"].append(evobs)
    galaxies[name]["Vgas"].append(vgas)
    galaxies[name]["Vdisk"].append(vdisk)
    galaxies[name]["Vbul"].append(vbul)

for g in galaxies:
    for k in galaxies[g]:
        if isinstance(galaxies[g][k], list):
            galaxies[g][k] = np.array(galaxies[g][k])

print(f"Cargadas {len(galaxies)} galaxias")

# -----------------------------------------------------------------------
# Modelos
# -----------------------------------------------------------------------
def make_exp_model(g):
    """v_obs^2 = Y*(v_disk^2+v_bul^2) + v_gas^2 + V_asy^2*(1-exp(-r/r_c))"""
    r = g["R"]; vgas = g["Vgas"]; vdisk = g["Vdisk"]; vbul = g["Vbul"]
    def model(r2, Y, V_asy, r_c):
        vbar2 = Y * (vdisk**2 + vbul**2) + vgas**2
        vvac2 = V_asy**2 * (1.0 - np.exp(-np.maximum(r2, 0.0) / np.maximum(r_c, 0.01)))
        return np.sqrt(np.maximum(vbar2 + vvac2, 0.1))
    return model

def make_laplace_model(g):
    """v_obs^2 = Y*(v_disk^2+v_bul^2) + v_gas^2 + V_asy^2 - M/r"""
    r = g["R"]; vgas = g["Vgas"]; vdisk = g["Vdisk"]; vbul = g["Vbul"]
    def model(r2, Y, V_asy, M):
        vbar2 = Y * (vdisk**2 + vbul**2) + vgas**2
        vvac2 = V_asy**2 - M / np.maximum(r2, 0.1)
        return np.sqrt(np.maximum(vbar2 + vvac2, 0.1))
    return model

def fit_galaxia(g, model_maker, n_params=3):
    r = g["R"]; vobs = g["Vobs"]; evobs = g["eVobs"]
    V0 = np.median(vobs[r > 0.5*r.max()]) if len(r) > 3 else 100.0
    Y0 = 0.5; rc0 = 2.0; M0 = V0**2 * 1.0
    p0 = [Y0, V0, rc0] if model_maker == make_exp_model else [Y0, V0, M0]
    model = model_maker(g)
    lb = [0.01, 0.1, 0.01] if model_maker == make_exp_model else [0.01, 0.1, -1e6]
    ub = [5.0, 500.0, 50.0] if model_maker == make_exp_model else [5.0, 500.0, 1e6]
    try:
        popt, pcov = curve_fit(model, r, vobs, p0=p0, sigma=evobs, maxfev=5000,
                               bounds=(lb, ub), ftol=1e-6, xtol=1e-6)
        perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(np.diag(pcov))) else [0]*n_params
        residuals = vobs - model(r, *popt)
        chi2 = np.sum((residuals / np.maximum(evobs, 0.1))**2)
        chi2_red = chi2 / max(len(r) - n_params, 1)
        return {"Y": popt[0], "V_asy": popt[1],
                "rc": popt[2] if model_maker == make_exp_model else 0,
                "M": popt[2] if model_maker == make_laplace_model else 0,
                "e_Y": perr[0], "e_V": perr[1], "e_3": perr[2],
                "chi2_red": chi2_red, "N": len(r)}, popt
    except Exception as e:
        return None, None

# -----------------------------------------------------------------------
# Fit exponencial a TODAS las 175 galaxias
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("AJUSTE EXPONENCIAL v_vac^2 = V_asy^2*(1 - exp(-r/r_c))")
print("Todas las galaxias SPARC")
print("=" * 75)
print(f"\n{'Galaxia':<14}{'Y':<10}{'V_asy':<10}{'r_c':<10}{'chi2_red':<10}{'N':<6}{'Status':<10}")
print("-" * 70)

resultados_exp = {}
t0 = time.time()
n_ok = 0; n_fail = 0
gal_names = sorted(galaxies.keys())
for gname in gal_names:
    g = galaxies[gname]
    if len(g["R"]) < 4:
        resultados_exp[gname] = None
        n_fail += 1
        continue
    res, _ = fit_galaxia(g, make_exp_model)
    resultados_exp[gname] = res
    if res and res["chi2_red"] < 100:
        n_ok += 1
        status = "OK" if res["chi2_red"] < 5 else "alta"
        print(f"{gname:<14}{res['Y']:<10.3f}{res['V_asy']:<10.1f}{res['rc']:<10.2f}{res['chi2_red']:<10.2f}{res['N']:<6}{status:<10}")
    else:
        n_fail += 1
        print(f"{gname:<14}{'---':<10}{'---':<10}{'---':<10}{'---':<10}{'---':<6}{'FAIL':<10}")

t_total = time.time() - t0
chi2s = [r["chi2_red"] for r in resultados_exp.values() if r is not None and r["chi2_red"] < 100]

print(f"\nRESUMEN EXPONENCIAL:")
print(f"  Galaxias OK: {n_ok}/{len(galaxies)}")
print(f"  Tiempo: {t_total:.1f}s")
print(f"  chi2_red medio (todas): {np.mean(chi2s):.2f} +/- {np.std(chi2s):.2f}")
print(f"  chi2_red mediana: {np.median(chi2s):.2f}")
print(f"  chi2_red < 1: {sum(1 for c in chi2s if c < 1)}")
print(f"  chi2_red < 2: {sum(1 for c in chi2s if c < 2)}")
print(f"  chi2_red < 3: {sum(1 for c in chi2s if c < 3)}")
print(f"  chi2_red < 5: {sum(1 for c in chi2s if c < 5)}")

# -----------------------------------------------------------------------
# Fit Laplace a TODAS las 175 galaxias (comparacion)
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("AJUSTE LAPLACE (comparacion) v_vac^2 = V_asy^2 - M/r")
print("=" * 75)

resultados_lap = {}
n_ok_l = 0; n_fail_l = 0
for gname in gal_names:
    g = galaxies[gname]
    if len(g["R"]) < 4:
        resultados_lap[gname] = None
        n_fail_l += 1
        continue
    res, _ = fit_galaxia(g, make_laplace_model)
    resultados_lap[gname] = res
    if res and res["chi2_red"] < 100:
        n_ok_l += 1
    else:
        n_fail_l += 1

chi2s_lap = [r["chi2_red"] for r in resultados_lap.values() if r is not None and r["chi2_red"] < 100]

print(f"\nRESUMEN LAPLACE:")
print(f"  Galaxias OK: {n_ok_l}/{len(galaxies)}")
print(f"  chi2_red medio: {np.mean(chi2s_lap):.2f} +/- {np.std(chi2s_lap):.2f}")
print(f"  chi2_red mediana: {np.median(chi2s_lap):.2f}")
print(f"  chi2_red < 1: {sum(1 for c in chi2s_lap if c < 1)}")
print(f"  chi2_red < 2: {sum(1 for c in chi2s_lap if c < 2)}")
print(f"  chi2_red < 3: {sum(1 for c in chi2s_lap if c < 3)}")
print(f"  chi2_red < 5: {sum(1 for c in chi2s_lap if c < 5)}")

# -----------------------------------------------------------------------
# Comparacion directa
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("COMPARACION DIRECTA EXP vs LAP (galaxias con ambos ajustes OK)")
print("=" * 75)
print(f"\n{'Galaxia':<14}{'Exp_chi2':<10}{'Lap_chi2':<10}{'mejor':<10}")
print("-" * 44)
mejor_exp = 0; mejor_lap = 0
for gname in gal_names:
    re = resultados_exp.get(gname)
    rl = resultados_lap.get(gname)
    if re and rl and re["chi2_red"] < 100 and rl["chi2_red"] < 100:
        mejor = "EXP" if re["chi2_red"] < rl["chi2_red"] else "LAP"
        if mejor == "EXP": mejor_exp += 1
        else: mejor_lap += 1
        print(f"{gname:<14}{re['chi2_red']:<10.2f}{rl['chi2_red']:<10.2f}{mejor:<10}")

print(f"\nEXP gana: {mejor_exp}, LAP gana: {mejor_lap}")

# -----------------------------------------------------------------------
# Top 10 mejores galaxias (EXP)
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("TOP 10 MEJORES AJUSTES EXPONENCIAL")
print("=" * 75)
sorted_exp = sorted([(g, r) for g, r in resultados_exp.items() if r and r["chi2_red"] < 100],
                    key=lambda x: x[1]["chi2_red"])
print(f"\n{'Galaxia':<14}{'Y':<10}{'V_asy':<10}{'r_c':<10}{'chi2_red':<10}{'N':<6}")
print("-" * 60)
for gname, res in sorted_exp[:10]:
    print(f"{gname:<14}{res['Y']:<10.3f}{res['V_asy']:<10.1f}{res['rc']:<10.2f}{res['chi2_red']:<10.2f}{res['N']:<6}")

# -----------------------------------------------------------------------
# NGC 3198 detalle
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("NGC 3198 DETALLE (EXP vs LAP vs OBS)")
print("=" * 75)
g = galaxies["NGC3198"]
r = g["R"]; vobs = g["Vobs"]; evobs = g["eVobs"]
vgas = g["Vgas"]; vdisk = g["Vdisk"]; vbul = g["Vbul"]

_, popt_exp = fit_galaxia(g, make_exp_model)
_, popt_lap = fit_galaxia(g, make_laplace_model)

if popt_exp is not None and popt_lap is not None:
    Ye, Ve, rc = popt_exp; Yl, Vl, Ml = popt_lap
    model_e = make_exp_model(g); model_l = make_laplace_model(g)
    vexp = model_e(r, *popt_exp); vlap = model_l(r, *popt_lap)
    vbar = np.sqrt(Ye * (vdisk**2 + vbul**2) + vgas**2)
    vvac = np.sqrt(np.maximum(Ve**2 * (1.0 - np.exp(-r/np.maximum(rc, 0.01))), 0.1))

    print(f"\n{'r':>5} {'Vobs':>6} {'eV':>5} {'Vbar':>6} {'Vvac':>6} {'Vexp':>6} {'Vlap':>6} {'O-C_E':>7} {'O-C_L':>7}")
    print("-" * 62)
    for i in range(len(r)):
        oce = vobs[i] - vexp[i]; ocl = vobs[i] - vlap[i]
        print(f"{r[i]:5.2f} {vobs[i]:6.1f} {evobs[i]:5.1f} {vbar[i]:6.1f} {vvac[i]:6.1f} {vexp[i]:6.1f} {vlap[i]:6.1f} {oce:7.1f} {ocl:7.1f}")

    re = resultados_exp.get("NGC3198")
    rl = resultados_lap.get("NGC3198")
    if re and rl:
        print(f"\nEXP: Y={Ye:.3f}, V_asy={Ve:.0f}, r_c={rc:.2f}, chi2={re['chi2_red']:.2f}")
        print(f"LAP: Y={Yl:.3f}, V_asy={Vl:.0f}, M={Ml:.0f}, chi2={rl['chi2_red']:.2f}")

# -----------------------------------------------------------------------
# Distribucion de V_asy y r_c
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("DISTRIBUCION DE PARAMETROS (EXPONENCIAL)")
print("=" * 75)
V_vals = [r["V_asy"] for r in resultados_exp.values() if r and r["chi2_red"] < 10]
rc_vals = [r["rc"] for r in resultados_exp.values() if r and r["chi2_red"] < 10]
Y_vals = [r["Y"] for r in resultados_exp.values() if r and r["chi2_red"] < 10]

if V_vals:
    print(f"\nV_asy: media={np.mean(V_vals):.0f} +/- {np.std(V_vals):.0f} km/s, mediana={np.median(V_vals):.0f}")
    print(f"r_c:   media={np.mean(rc_vals):.2f} +/- {np.std(rc_vals):.2f} kpc, mediana={np.median(rc_vals):.2f}")
    print(f"Y:     media={np.mean(Y_vals):.3f} +/- {np.std(Y_vals):.3f}, mediana={np.median(Y_vals):.3f}")

# -----------------------------------------------------------------------
# Conclusion
# -----------------------------------------------------------------------
print("\n" + "=" * 75)
print("CONCLUSION")
print("=" * 75)
print(f"""
  Perfil exponencial resuelve el problema de signo del perfil de Laplace:
  v_vac^2(r) = V_asy^2 * (1 - exp(-r/r_c))

  Ventajas:
  1. v_vac^2(0) = 0 (no contribucion en el centro)
  2. v_vac^2(r) -> V_asy^2 (asintotico constante, consistente con curvas planas)
  3. SIEMPRE positiva (no hay repulsion)
  4. Un solo parametro adicional: r_c (escala de transicion)

  Resultados ({n_ok} galaxias):
  - chi2_red medio = {np.mean(chi2s):.2f}
  - Galaxias con chi2_red < 2: {sum(1 for c in chi2s if c < 2)}/{n_ok}
  - Galaxias con chi2_red < 3: {sum(1 for c in chi2s if c < 3)}/{n_ok}
  - Comparacion: EXP gana a LAP en {mejor_exp} vs {mejor_lap} galaxias

  Proximo paso: comparar con modelo de materia oscura NFW
  para las 175 galaxias y verificar que chi2_exp ~ chi2_nfw.
""")
