"""
Avbay_Steady.py
1D steady heat balance of the sealed avbay: PCB temperatures and ANSYS boundary conditions
for worst, median and low pad conditions.
"""

import numpy as np
from scipy.optimize import brentq, fsolve

SIG = 5.670e-8                    # Stefan-Boltzmann constant [W/m^2 K^4]

# PCBs: one-face area (KiCad) [m^2], Ansys body volume (Details > Volume) [m^3]
NAMES = ["C6", "G6", "ENIAC"]
A = 2 * np.array([0.01451, 0.01451, 0.00587])          # both faces [m^2]
V = np.array([2.1857e-5, 2.1876e-5, 9.3937e-6])        # [m^3]

BETA_DEG = 45.0                   # sun elevation, assumed [deg]
RHO_G = 0.40                      # ground albedo, desert sand (Wikipedia, Albedo) [-]
H_I = 3.0                         # inside convection, assumed; Churchill-Chu gives 2.1-5.5 [W/m^2 K]
A_O, A_W = 0.245, 0.268           # outer tube, inner wall area [m^2]
L_WALL = 0.175 * 0.0254           # wall thickness [m]
K_FR4 = 0.29                      # FR4 through-plane conductivity (Wikipedia, FR-4) [W/m K]
A_PAINT, E_PAINT = 0.60, 0.90     # paint absorptivity, emissivity, assumed until paint is chosen [-]
E_FR4 = 0.90                      # FR4 emissivity, non-metal (FLIR) [-]


def C2K(T): return T + 273.15
def K2C(T): return T - 273.15
def CF(T): return f"{K2C(T):.1f} C ({K2C(T) * 1.8 + 32:.1f} F)"


def solve(T_amb_c, T_gnd_c, G_b, G_d, h_o, P):
    """Skin, wall, air and PCB temperatures [K] for one set of pad conditions."""
    P = np.array(P)
    T_inf = C2K(T_amb_c)
    T_sky = 0.0552 * T_inf**1.5                            # Swinbank
    T_rad = ((T_sky**4 + C2K(T_gnd_c)**4) / 2) ** 0.25     # half sky, half ground
    b = np.radians(BETA_DEG)
    G_h = G_b * np.sin(b) + G_d                            # global horizontal
    # absorbed sun: beam on cylinder + isotropic sky + ground-reflected
    q_sun = A_PAINT * (G_b * np.cos(b) / np.pi + G_d / 2 + RHO_G * G_h / 2)

    def q_out(Ts):
        """Newton cooling + Stefan-Boltzmann to T_rad [W/m^2]."""
        return h_o * (Ts - T_inf) + E_PAINT * SIG * (Ts**4 - T_rad**4)

    def air(Tp, Tw):
        """Sealed-air balance [K]."""
        return (A @ Tp + A_W * Tw) / (A.sum() + A_W)

    def board_loss(Tp, Tw):
        """Convection to air + gray-body radiation via wall radiosity [W]."""
        w = E_FR4 * A / A_W
        J_w = SIG * (E_FR4 * Tw**4 + (1 - E_FR4) * (w @ Tp**4)) / (E_FR4 + (1 - E_FR4) * w.sum())
        return H_I * A * (Tp - air(Tp, Tw)) + E_FR4 * A * (SIG * Tp**4 - J_w)

    # skin: sun + PCB heat = convection + radiation out
    T_s = brentq(lambda T: q_sun * A_O + P.sum() - A_O * q_out(T), T_inf - 50, T_inf + 300)
    T_w = T_s + P.sum() * L_WALL / (K_FR4 * A_W)           # Fourier's law
    T_p, _, ok, msg = fsolve(lambda Tp: P - board_loss(Tp, T_w), np.full(len(P), T_w + 20),
                             xtol=1e-12, full_output=True)
    if ok != 1:
        raise RuntimeError(f"PCB solve failed: {msg}")
    return dict(P=P, h_o=h_o, T_inf=T_inf, T_rad=T_rad, q_sun=q_sun,
                T_s=T_s, T_w=T_w, T_a=air(T_p, T_w), T_p=T_p)


def report(title, r):
    """Temperatures and the matching ANSYS boundary conditions."""
    Q_GEN = r["P"] / V                                     # q''' = P / V_body [W/m^3]
    print(f"\n{title}")
    print(f"Ambient: {CF(r['T_inf'])}")
    print(f"Internal Air: {CF(r['T_a'])}")
    print("PCB      Heat [W]   IHG [W/m^3]   Temp")
    for n, p, q, T in zip(NAMES, r["P"], Q_GEN, r["T_p"]):
        print(f"{n:<8}{p:8.2f}{q:14,.0f}   {CF(T)}")
    print("ANSYS STEADY-STATE BOUNDARY CONDITIONS")
    bcs = [
        ("Initial Temperature", "all bodies", CF(r["T_inf"])),
        ("Heat Flux", "outer shell", f"{r['q_sun']:.1f} W/m^2"),
        ("Convection", "outer shell", f"h {r['h_o']:g} W/m^2C, ambient {CF(r['T_inf'])}"),
        ("Radiation", "outer shell", f"To Ambient, e {E_PAINT:g}, ambient {CF(r['T_rad'])}"),
        *[("Internal Heat Gen" if i == 0 else "", n, f"{q:,.0f} W/m^3")
          for i, (n, q) in enumerate(zip(NAMES, Q_GEN))],
        ("Convection", "inside faces", f"h {H_I:g} W/m^2C, ambient {CF(r['T_a'])}"),
        ("Radiation", "inside faces", f"Surface to Surface, Perfect, Enclosure 1, ambient {CF(r['T_w'])}"),
        ("", "", f"e: FR4 {E_FR4:g}"),
    ]
    for bc, scope, value in bcs:
        print(f"  {bc:<20}{scope:<14}{value}")


# Worst case: DTEG max air, 3 mph wind, full sun, ENIAC transmitting at 8.4 V
T_AMB_C = 42.0                    # air, DTEG App. G max (108 F) [C]
T_GND_C = 60.0                    # ground, placeholder until NASA launch-site data [C]
G_B, G_D = 900.0, 100.0           # direct (ASTM G173 AM1.5 direct), diffuse (assumed) sun [W/m^2]
H_O = 10.0                        # outside convection, assumed 3 mph wind (Churchill-Bernstein) [W/m^2 K]
P_PCB = [4.72649, 0.69056, 6.64315]   # C6, G6, ENIAC heat [W]
report("WORST CASE", solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, P_PCB))

# Median: DTEG average high, average wind, full sun, ENIAC idle at nominal 7.4 V
T_AMB_C = 34.0                    # air, DTEG App. G average high 33-35 C [C]
T_GND_C = 52.0                    # ground, same 18 C over air as worst case, placeholder [C]
G_B, G_D = 900.0, 100.0           # direct (ASTM G173 AM1.5 direct), diffuse (assumed) sun [W/m^2]
H_O = 21.6                        # outside convection, 11.6 mph DTEG average (Churchill-Bernstein) [W/m^2 K]
P_PCB = [4.72649, 0.69056, 7.4 * 0.65767]   # ENIAC idle: 7.4 V x 0.658 A [W]
report("MEDIAN CASE", solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, P_PCB))

# Low: DTEG average low, no sun, average wind, ENIAC idle at nominal 7.4 V
T_AMB_C = 20.5                    # air, DTEG App. G average low 19-22 C [C]
T_GND_C = 20.5                    # ground, no solar heating, taken at air temperature [C]
G_B, G_D = 0.0, 0.0               # no sun [W/m^2]
H_O = 22.0                        # outside convection, 11.6 mph DTEG average (Churchill-Bernstein) [W/m^2 K]
P_PCB = [4.72649, 0.69056, 7.4 * 0.65767]   # ENIAC idle: 7.4 V x 0.658 A [W]
report("LOW CASE", solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, P_PCB))
