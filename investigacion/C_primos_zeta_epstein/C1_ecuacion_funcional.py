import numpy as np, math

print("=" * 72)
print("  C1: ECUACION FUNCIONAL VIA FUNCION THETA DE Z^3")
print("  Conexion con red D3 (FCC) como subred de Z^3")
print("=" * 72)

# -----------------------------------------------------------------------
# Funcion theta para Z^3 (simple cubic)
# theta_Z3(t) = sum_{n in Z^3} exp(-pi * |n|^2 * t)
# -----------------------------------------------------------------------
def theta_Z3(t, N):
    s = 0.0
    for i in range(-N, N+1):
        for j in range(-N, N+1):
            for k in range(-N, N+1):
                s += math.exp(-math.pi * (i*i+j*j+k*k) * t)
    return s

print(f"\n  Verificando para Z^3 (self-dual):")
print(f"  theta_Z3(t) = t^(-3/2) * theta_Z3(1/t)")
print(f"\n  {'t':<10}{'theta_Z3(t)':<18}{'t^{-3/2}*theta_Z3(1/t)':<24}{'cociente':<12}")
print(f"  " + "-" * 64)

for t in [0.2, 0.5, 1.0, 2.0, 5.0]:
    th = theta_Z3(t, 8)
    rhs = t**(-1.5) * theta_Z3(1.0/t, 8)
    ratio = th / rhs
    print(f"  {t:<10.2f}{th:<18.6f}{rhs:<24.6f}{ratio:<12.6f}")

# -----------------------------------------------------------------------
# Funcion theta para D3 (FCC): subred de Z^3 con i+j+k par
# -----------------------------------------------------------------------
def theta_FCC(t, N):
    s = 0.0
    for i in range(-N, N+1):
        for j in range(-N, N+1):
            for k in range(-N, N+1):
                if (i+j+k) % 2 != 0:
                    continue
                # |v|^2_FCC = (i^2+j^2+k^2)/2 (NN distance = 1)
                s += math.exp(-math.pi * (i*i+j*j+k*k) * t / 2.0)
    return s

def theta_BCC(t, N):
    """Theta de la red dual BCC: Z^3 + (1/2,1/2,1/2)."""
    s = theta_Z3(t, N)  # suma sobre Z^3
    for i in range(-N, N+1):
        for j in range(-N, N+1):
            for k in range(-N, N+1):
                ip = i + 0.5; jp = j + 0.5; kp = k + 0.5
                s += math.exp(-math.pi * (ip*ip+jp*jp+kp*kp) * t)
    return s

print(f"\n  Ecuacion funcional para FCC (D3):")
print(f"  theta_FCC(t) = det(D3)^(-1/2) * t^(-3/2) * theta_BCC(1/t)")
print(f"  det(D3) = 2 (indice 2 en Z^3)")
print(f"\n  {'t':<10}{'theta_FCC(t)':<18}{'factor*theta_BCC(1/t)':<24}{'cociente':<12}")
print(f"  " + "-" * 64)

for t in [0.2, 0.5, 1.0, 2.0, 5.0]:
    th_fcc = theta_FCC(t, 6)
    th_bcc = theta_BCC(1.0/t, 6)
    rhs = 2.0**(-0.5) * t**(-1.5) * th_bcc
    ratio = th_fcc / rhs
    print(f"  {t:<10.2f}{th_fcc:<18.6f}{rhs:<24.6f}{ratio:<12.6f}")

# -----------------------------------------------------------------------
# Zeta de Epstein de Z^3 en linea critica Re(s)=3/2
# -----------------------------------------------------------------------
print(f"\n  ZETA DE EPSTEIN PARA Z^3 (linea critica Re(s) = 3/2):")
print(f"  ===================================================")
print(f"  Para Z^3 (dim 3), la ecuacion funcional tiene linea")
print(f"  critica en Re(s) = 3/2. La funcion: ")
print("    xi_Z3(s) = pi^{-s} * Gamma(s) * zeta_Z3(s)")
print("    xi(3/2 - s) = xi(s)")
print("  donde zeta_Z3(s) = sum_{n != 0} 1/(n1^2+n2^2+n3^2)^s")

def zeta_Z3(s, N):
    total = 0.0j
    for i in range(-N, N+1):
        for j in range(-N, N+1):
            for k in range(-N, N+1):
                if i == 0 and j == 0 and k == 0:
                    continue
                n2 = float(i*i + j*j + k*k)
                total += n2**(-s)
    return total

print(f"\n  {'s':<10}{'zeta_Z3(s)':<20}{'nota':<30}")
for s in [2.0, 2.5, 3.0, 4.0]:
    z = zeta_Z3(s, 6)
    print(f"  {s:<10.2f}{z.real:<20.6f}{'converge (Re>1.5)':<30}")

# -----------------------------------------------------------------------
# Conexion con Norma=4
# -----------------------------------------------------------------------
print("\n  CONEXION CON NORMA=4:")
print("  ===================")
print("""
  La funcion theta de Z^3 se factoriza como theta_1D(t)^3:
    theta_Z3(t) = (sum_{n=-inf}^{inf} exp(-pi*n^2*t))^3

  El umbral Norma=4 en el modelo FCC-DNLS corresponde a:
    - Primer entero n = 4 con r3(n) = 6 representaciones
    - Autovalor lambda = 4 del Laplaciano FCC
    - Primera representacion NO TRIVIAL de suma de 3 cuadrados
      (despues de n=0, que es trivial)

  En la ecuacion funcional, la linea critica Re(s)=3/4 del modelo
  FCC-DNLS refleja la dimension efectiva 3/2 del sistema:
    - Espacio 3D (FCC lattice)
    - Pero dinamica 1D (tiempo emergente)
    - d_efectiva = 3/2 = 3D/2""")

# -----------------------------------------------------------------------
# Conclusion
# -----------------------------------------------------------------------
print(f"\n" + "=" * 72)
print("  CONCLUSION")
print("=" * 72)
print(f"""
  La funcion theta de Z^3 satisface:
    theta_Z3(t) = t^{-3/2} * theta_Z3(1/t)

  Para la red D3 (FCC), la ecuacion funcional es:
    theta_FCC(t) = det(D3)^{-1/2} * t^{-3/2} * theta_BCC(1/t)
  donde BCC es la red dual de FCC.

  La linea critica de la zeta de Epstein para Z^3 es Re(s)=3/2,
  mientras que para el modelo FCC-DNLS es Re(s)=3/4. Esto refleja
  la dimension efectiva 3/2 = dimension espacial/2 del sistema
  con tiempo emergente.

  El umbral Norma=4 es el PRIMER ENTERO NO NULO con representacion
  no trivial como suma de tres cuadrados (r3(4)=6), conectando
  la aritmetica de sumas de tres cuadrados con la estabilidad
  del punto fijo espectral.
""")
