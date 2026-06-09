import math, sys
import numpy as np

# ------------------------------------------------------------
# Cargar datos SPARC
# ------------------------------------------------------------
data_raw = open(r"C:\Users\famil\Desktop\serie_armonica\SPARC\table2.dat", "rb").read().decode("latin-1")
lines = [l for l in data_raw.split("\n") if l.strip() and not l.startswith('\r')]

galaxies = {}
for line in lines:
    name = line[0:11].strip()
    try:
        dist = float(line[12:18].strip())
        rad = float(line[19:25].strip())
        vobs = float(line[26:32].strip())
        evobs = float(line[33:38].strip())
        vgas = float(line[39:45].strip())
        vdisk = float(line[46:52].strip())
        vbul = float(line[53:59].strip())
    except (ValueError, IndexError):
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
    galaxies[g]["R"] = np.array(galaxies[g]["R"])
    galaxies[g]["Vobs"] = np.array(galaxies[g]["Vobs"])
    galaxies[g]["eVobs"] = np.array(galaxies[g]["eVobs"])
    galaxies[g]["Vgas"] = np.array(galaxies[g]["Vgas"])
    galaxies[g]["Vdisk"] = np.array(galaxies[g]["Vdisk"])
    galaxies[g]["Vbul"] = np.array(galaxies[g]["Vbul"])

print(f"Galaxias cargadas: {len(galaxies)}")
print()

# ------------------------------------------------------------
# Modelo: v_model^2 = Y * (v_disk^2 + v_bul^2) + v_gas^2 + V_asy^2 - M/r
# ------------------------------------------------------------
def modelo_laplace(r, Y, V_asy, M):
    """Laplace vacuum tension model.
    v_obs^2 = Y*(v_disk^2+v_bul^2) + v_gas^2 + V_asy^2 - M/r
    """
    vbar2 = Y * (v_disk**2 + v_bul**2) + v_gas**2
    vvac2 = V_asy**2 - M / np.maximum(r, 0.01)
    vtot2 = vbar2 + vvac2
    return np.sqrt(np.maximum(vtot2, 0))

# ------------------------------------------------------------
# Ajuste para cada galaxia
# ------------------------------------------------------------
def fit_galaxia(gname, Y0=0.5, V0=None, M0=None):
    g = galaxies[gname]
    r = g["R"]
    vobs = g["Vobs"]
    evobs = g["eVobs"]
    vgas = g["Vgas"]
    vdisk = g["Vdisk"]
    vbul = g["Vbul"]
    w = 1.0 / np.maximum(evobs, 0.1)

    if V0 is None:
        V0 = np.median(vobs[r > 0.5*np.max(r)])
    if M0 is None:
        M0 = V0**2 * 0.5

    def chi2(params):
        Y, V, M = params
        vbar2 = Y * (vdisk**2 + vbul**2) + vgas**2
        vvac2 = V**2 - M / np.maximum(r, 0.01)
        vtot2 = vbar2 + vvac2
        vmod = np.sqrt(np.maximum(vtot2, 0))
        return np.sum(w * (vobs - vmod)**2)

    # Simple grid search
    best = (Y0, V0, M0)
    best_c2 = chi2(best)

    for Y in np.linspace(0.1, 2.0, 20):
        for V in np.linspace(0.5*V0, 1.5*V0, 20):
            for M in np.linspace(0, 2*M0, 20):
                c2 = chi2((Y, V, M))
                if c2 < best_c2:
                    best_c2 = c2
                    best = (Y, V, M)

    Yb, Vb, Mb = best
    N = len(r)
    chi2_red = best_c2 / max(N - 3, 1)
    r_transition = Mb / Vb**2 if Vb > 0 else 0

    return {
        "Y": Yb, "V_asy": Vb, "M": Mb,
        "chi2_red": chi2_red, "r_transition": r_transition,
        "N_points": N
    }

# ------------------------------------------------------------
# Galaxias a analizar
# ------------------------------------------------------------
targets = ["NGC3198", "NGC2403", "NGC2903", "NGC5055", "NGC6503", "NGC6946"]

print(f"{'Galaxia':<12}{'Y':<10}{'V_asy':<10}{'M':<12}{'r_trans':<10}{'chi2_red':<10}{'N':<6}")
print("-" * 70)
for gname in targets:
    if gname in galaxies:
        res = fit_galaxia(gname)
        print(f"{gname:<12}{res['Y']:<10.3f}{res['V_asy']:<10.1f}{res['M']:<12.1f}{res['r_transition']:<10.2f}{res['chi2_red']:<10.4f}{res['N_points']:<6}")
    else:
        print(f"{gname:<12}NO ENCONTRADA")

# ------------------------------------------------------------
# Detalle NGC 3198
# ------------------------------------------------------------
print("\n" + "=" * 60)
print("DETALLE NGC 3198")
print("=" * 60)
g = galaxies["NGC3198"]
r = g["R"]
vobs = g["Vobs"]
evobs = g["eVobs"]
vgas = g["Vgas"]
vdisk = g["Vdisk"]
vbul = g["Vbul"]
res = fit_galaxia("NGC3198")
Yb, Vb, Mb = res["Y"], res["V_asy"], res["M"]

vbar = np.sqrt(Yb * (vdisk**2 + vbul**2) + vgas**2)
vvac = np.sqrt(np.maximum(Vb**2 - Mb / np.maximum(r, 0.01), 0))
vmod = np.sqrt(np.maximum(Yb * (vdisk**2 + vbul**2) + vgas**2 + Vb**2 - Mb / np.maximum(r, 0.01), 0))

print(f"{'r(kpc)':<10}{'Vobs':<10}{'eVobs':<10}{'Vbar':<10}{'Vvac':<10}{'Vmod':<10}{'residuo':<10}")
print("-" * 70)
for i in range(len(r)):
    resid = vobs[i] - vmod[i]
    print(f"{r[i]:<10.2f}{vobs[i]:<10.1f}{evobs[i]:<10.1f}{vbar[i]:<10.1f}{vvac[i]:<10.1f}{vmod[i]:<10.1f}{resid:<10.1f}")

print(f"\nRESULTADOS NGC 3198:")
print(f"  Y (M/L) = {Yb:.3f}")
print(f"  V_asy = {Vb:.1f} km/s")
print(f"  M = {Mb:.1f} km^2/s^2 kpc")
print(f"  r_transition = {res['r_transition']:.2f} kpc")
print(f"  chi2_red = {res['chi2_red']:.4f} (N={res['N_points']})")
