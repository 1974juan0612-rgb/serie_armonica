"""
SPARC rotation curve fitting with Laplace vacuum profile.
Model: v_obs^2 = Y * (v_disk^2 + v_bul^2) + v_gas^2 + V_asy^2 - M_k / r
"""
import numpy as np
from scipy.optimize import curve_fit

# ------------------------------------------------------------
# Load SPARC data
# ------------------------------------------------------------
raw = open(r"C:\Users\famil\Desktop\serie_armonica\SPARC\table2.dat", "rb").read().decode("latin-1")
lines = [l for l in raw.split("\n") if l.strip() and len(l) >= 76]

galaxies = {}
for line in lines:
    try:
        name = line[0:11].strip()
        rad = float(line[19:25].strip())
        vobs = float(line[26:32].strip())
        evobs = float(line[33:38].strip())
        vgas = float(line[39:45].strip())
        vdisk = float(line[46:52].strip())
        vbul = float(line[53:59].strip())
    except ValueError:
        continue
    if name not in galaxies:
        galaxies[name] = {"R": [], "Vobs": [], "eVobs": [], "Vgas": [], "Vdisk": [], "Vbul": []}
    galaxies[name]["R"].append(rad)
    galaxies[name]["Vobs"].append(vobs)
    galaxies[name]["eVobs"].append(evobs)
    galaxies[name]["Vgas"].append(vgas)
    galaxies[name]["Vdisk"].append(vdisk)
    galaxies[name]["Vbul"].append(vbul)

for g in galaxies:
    for k in galaxies[g]:
        galaxies[g][k] = np.array(galaxies[g][k])

print(f"Loaded {len(galaxies)} galaxies\n")

# ------------------------------------------------------------
# Fitting with closures to avoid globals
# ------------------------------------------------------------
def make_laplace_model(g):
    r = g["R"]; vgas = g["Vgas"]; vdisk = g["Vdisk"]; vbul = g["Vbul"]
    def model(r2, Y, V_asy, M):
        vbar2 = Y * (vdisk**2 + vbul**2) + vgas**2
        vvac2 = V_asy**2 - M / np.maximum(r2, 0.1)
        return np.sqrt(np.maximum(vbar2 + vvac2, 0.1))
    return model

def make_nfw_model(g):
    r = g["R"]; vgas = g["Vgas"]; vdisk = g["Vdisk"]; vbul = g["Vbul"]
    def model(r2, Y, V_200, r_s):
        x = r2 / np.maximum(r_s, 0.1)
        vbar2 = Y * (vdisk**2 + vbul**2) + vgas**2
        vhalo2 = V_200**2 * (np.log(1+x) - x/(1+x)) / np.maximum(x, 0.01)
        return np.sqrt(vbar2 + vhalo2)
    return model

def fit_galaxy(gname, model_maker=make_laplace_model):
    g = galaxies[gname]
    r = g["R"]; vobs = g["Vobs"]; evobs = g["eVobs"]

    V0 = np.median(vobs[r > 0.5*r.max()])
    M0 = V0**2 * 1.0
    Y0 = 0.5

    model = model_maker(g)
    try:
        popt, pcov = curve_fit(model, r, vobs, p0=[Y0, V0, M0], sigma=evobs, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        residuals = vobs - model(r, *popt)
        chi2 = np.sum((residuals / evobs)**2)
        chi2_red = chi2 / max(len(r) - 3, 1)
        r_trans = popt[2] / popt[1]**2 if popt[1] > 0 else 0
        return {"Y": popt[0], "V_asy": popt[1], "M": popt[2],
                "e_Y": perr[0], "e_V": perr[1], "e_M": perr[2],
                "r_trans": r_trans, "chi2_red": chi2_red, "N": len(r)}, popt
    except Exception as e:
        print(f"  {gname} FAILED: {e}")
        return None, None

# ------------------------------------------------------------
# Fit galaxies
# ------------------------------------------------------------
targets = ["NGC3198", "NGC2403", "NGC2903", "NGC5055", "NGC6503", "NGC6946"]

print("LAPLACE PROFILE VACUUM TENSION")
print(f"{'Galaxy':<12}{'Y':<10}{'V_asy':<10}{'M_corr':<12}{'r_trans':<10}{'chi2_red':<10}{'N':<6}")
print("-" * 70)
for gname in targets:
    if gname in galaxies:
        res, _ = fit_galaxy(gname, make_laplace_model)
        if res:
            print(f"{gname:<12}{res['Y']:<10.3f}{res['V_asy']:<10.1f}{res['M']:<12.0f}{res['r_trans']:<10.2f}{res['chi2_red']:<10.2f}{res['N']:<6}")

print("\nNFW DARK MATTER HALO (comparison)")
print(f"{'Galaxy':<12}{'Y':<10}{'V_200':<10}{'r_s':<10}{'chi2_red':<10}{'N':<6}")
print("-" * 58)
for gname in targets:
    if gname in galaxies:
        g = galaxies[gname]
        V200_guess = np.median(g["Vobs"][g["R"] > 0.5*g["R"].max()])
        res, _ = fit_galaxy(gname, make_nfw_model)
        if res:
            print(f"{gname:<12}{res['Y']:<10.3f}{res['V_asy']:<10.1f}{res['M']:<10.1f}{res['chi2_red']:<10.2f}{res['N']:<6}")

# ------------------------------------------------------------
# Detailed NGC 3198 table
# ------------------------------------------------------------
print("\n" + "=" * 75)
print("NGC 3198 - DETAIL")
print("=" * 75)
g = galaxies["NGC3198"]
r, vobs, evobs, vgas, vdisk, vbul = [g[k] for k in ["R","Vobs","eVobs","Vgas","Vdisk","Vbul"]]

res_lap, popt_lap = fit_galaxy("NGC3198", make_laplace_model)
res_nfw, popt_nfw = fit_galaxy("NGC3198", make_nfw_model)

Yl, Vl, Ml = popt_lap
Yn, Vn, rsn = popt_nfw
model_lap = make_laplace_model(g)
model_nfw = make_nfw_model(g)
vmod_lap = model_lap(r, *popt_lap)
vmod_nfw = model_nfw(r, *popt_nfw)

vbar = np.sqrt(Yl * (vdisk**2 + vbul**2) + vgas**2)
vvac = np.sqrt(np.maximum(Vl**2 - Ml/np.maximum(r,0.1), 0.1))

print(f"\n{'r':>5} {'Vobs':>6} {'eV':>5} {'Vbar':>6} {'Vvac':>6} {'Vlap':>6} {'Vnfw':>6} {'O-C_L':>7} {'O-C_N':>7}")
print("-" * 60)
for i in range(len(r)):
    ocl = vobs[i] - vmod_lap[i]
    ocn = vobs[i] - vmod_nfw[i]
    print(f"{r[i]:5.1f} {vobs[i]:6.0f} {evobs[i]:5.0f} {vbar[i]:6.0f} {vvac[i]:6.0f} {vmod_lap[i]:6.0f} {vmod_nfw[i]:6.0f} {ocl:7.1f} {ocn:7.1f}")

print(f"\nLAPLACE: Y={Yl:.3f}+-{res_lap['e_Y']:.3f}, V_asy={Vl:.0f}+-{res_lap['e_V']:.0f} km/s, M={Ml:.0f}+-{res_lap['e_M']:.0f}, r_t={res_lap['r_trans']:.2f} kpc, chi2={res_lap['chi2_red']:.2f}")
print(f"NFW:    Y={Yn:.3f}+-{res_nfw['e_Y']:.3f}, V_200={Vn:.0f}+-{res_nfw['e_V']:.0f} km/s, r_s={rsn:.1f}+-{res_nfw['e_M']:.1f} kpc, chi2={res_nfw['chi2_red']:.2f}")
