"""Pre-launch pad-soak cases: electronics OFF, ground heating ON."""
import numpy as np

from Avbay_Steady import NAMES, V, CF, Settings, solve

ELECTRONICS_OFF_W = np.zeros(3)


def report(title, r):
    """Pad temperatures and matching ANSYS boundary conditions."""
    cfg = r['settings']
    Q_GEN = r['P'] / V

    print(f'\n{title}')
    print(f"Ambient: {CF(r['T_inf'])}")
    print(f"Internal Air: {CF(r['T_a'])}")
    print('PCB      Heat [W]   IHG [W/m^3]   Temp')
    for n, p, q, T in zip(NAMES, r['P'], Q_GEN, r['T_p']):
        print(f'{n:<8}{p:8.2f}{q:14,.0f}   {CF(T)}')

    print('ANSYS STEADY-STATE BOUNDARY CONDITIONS')
    bcs = [
        ('Initial Temperature', 'all bodies', CF(r['T_inf'])),
        ('Heat Flux', 'outer shell', f"{r['q_sun']:.1f} W/m^2"),
        ('Convection', 'outer shell', f"h {r['h_o']:g} W/m^2C, ambient {CF(r['T_inf'])}"),
        ('Radiation', 'outer shell', f"e {cfg.exterior_emissivity:g}, ambient {CF(r['T_rad'])}"),
        *[('Internal Heat Gen' if i == 0 else '', n, f'{q:,.0f} W/m^3')
          for i, (n, q) in enumerate(zip(NAMES, Q_GEN))],
        ('Convection', 'inside faces', f"h {cfg.h_board:g} W/m^2C, ambient {CF(r['T_a'])}"),
        ('Radiation', 'inside faces', f"Surface to Surface, e {cfg.board_emissivity:g}, wall {CF(r['T_w'])}"),
    ]
    for bc, scope, value in bcs:
        print(f'  {bc:<20}{scope:<14}{value}')


# Hot pad: electronics OFF
T_AMB_C = 42.0
T_GND_C = 60.0
G_B, G_D = 900.0, 100.0
H_O = 10.0
report('HOT PAD - ELECTRONICS OFF',
       solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, ELECTRONICS_OFF_W))

# Moderate pad: electronics OFF
T_AMB_C = 34.0
T_GND_C = 52.0
G_B, G_D = 900.0, 100.0
H_O = 21.6
report('MODERATE PAD - ELECTRONICS OFF',
       solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, ELECTRONICS_OFF_W))

# Cool pad / no sun: electronics OFF
T_AMB_C = 20.5
T_GND_C = 20.5
G_B, G_D = 0.0, 0.0
H_O = 22.0
report('COOL / NO SUN - ELECTRONICS OFF',
       solve(T_AMB_C, T_GND_C, G_B, G_D, H_O, ELECTRONICS_OFF_W))
