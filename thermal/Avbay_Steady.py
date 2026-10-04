"""
Avbay_Steady.py
1D steady heat balance of the sealed avbay: PCB temperatures and ANSYS boundary conditions.
"""

import numpy as np
from scipy.optimize import brentq, fsolve

SIG = 5.670e-8                    # Stefan-Boltzmann constant [W/m^2 K^4]

# PCBs: heat [W], one-face area (KiCad) [m^2], Ansys body volume (Details > Volume) [m^3]
PCB = {"C6": (4.57134, 0.01451, 2.1857e-5),
       "G6": (0.53077, 0.01451, 2.1876e-5),
       "ENIAC": (6.8786, 0.00587, 9.3937e-6)}

T_AMB_C = 42.0                    # air [C]
T_GND_C = 60.0                    # ground [C]
G_B, G_D = 1000.0, 100.0          # direct, diffuse sun [W/m^2]
BETA_DEG = 45.0                   # sun elevation [deg]
RHO_G = 0.30                      # ground albedo [-]
H_O, H_I = 10.0, 3.0              # outside, inside convection [W/m^2 K]
A_O, A_W = 0.245, 0.268           # outer tube, inner wall area [m^2]
L_WALL = 0.175 * 0.0254           # wall thickness [m]
K_FR4 = 0.29                      # FR4 conductivity [W/m K]
A_PAINT, E_PAINT = 0.60, 0.90     # paint absorptivity, emissivity [-]
E_FR4, E_AL, E_STEEL = 0.90, 0.20, 0.35   # emissivity [-]


def C2K(T): return T + 273.15
def K2C(T): return T - 273.15
def CF(T): return f"{K2C(T):.1f} C ({K2C(T) * 1.8 + 32:.1f} F)"


NAMES = list(PCB)
P = np.array([PCB[n][0] for n in NAMES])               # [W]
A = np.array([2 * PCB[n][1] for n in NAMES])           # both faces [m^2]
V = np.array([PCB[n][2] for n in NAMES])               # [m^3]
Q_GEN = P / V                                          # q''' = P / V_body [W/m^3]

T_inf = C2K(T_AMB_C)
T_sky = 0.0552 * T_inf**1.5                            # Swinbank
T_rad = ((T_sky**4 + C2K(T_GND_C)**4) / 2) ** 0.25     # half sky, half ground
b = np.radians(BETA_DEG)
G_h = G_B * np.sin(b) + G_D                            # global horizontal
# absorbed sun: beam on cylinder + isotropic sky + ground-reflected
q_sun = A_PAINT * (G_B * np.cos(b) / np.pi + G_D / 2 + RHO_G * G_h / 2)


def q_out(Ts):
    """Newton cooling + Stefan-Boltzmann to T_rad [W/m^2]."""
    return H_O * (Ts - T_inf) + E_PAINT * SIG * (Ts**4 - T_rad**4)


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
T_a = air(T_p, T_w)

print(f"Ambient: {CF(T_inf)}")
print(f"Internal Air: {CF(T_a)}")
print("PCB      Heat [W]   IHG [W/m^3]   Temp")
for n, p, q, T in zip(NAMES, P, Q_GEN, T_p):
    print(f"{n:<8}{p:8.2f}{q:14,.0f}   {CF(T)}")

print("\nANSYS STEADY-STATE BOUNDARY CONDITIONS")
bcs = [
    ("Initial Temperature", "all bodies", CF(T_inf)),
    ("Heat Flux", "outer shell", f"{q_sun:.1f} W/m^2"),
    ("Convection", "outer shell", f"h {H_O:g} W/m^2C, ambient {CF(T_inf)}"),
    ("Radiation", "outer shell", f"To Ambient, e {E_PAINT:g}, ambient {CF(T_rad)}"),
    *[("Internal Heat Gen" if i == 0 else "", n, f"{q:,.0f} W/m^3")
      for i, (n, q) in enumerate(zip(NAMES, Q_GEN))],
    ("Convection", "inside faces", f"h {H_I:g} W/m^2C, ambient {CF(T_a)}"),
    ("Radiation", "inside faces", "Surface to Surface, Perfect, Enclosure 1"),
    ("", "", f"e: FR4 {E_FR4:g}, Al {E_AL:g}, steel {E_STEEL:g}"),
]
for bc, scope, value in bcs:
    print(f"  {bc:<20}{scope:<14}{value}")
