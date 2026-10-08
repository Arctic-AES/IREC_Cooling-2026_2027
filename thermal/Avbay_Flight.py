"""In-flight steady thermal sensitivity cases: electronics ON, ground heating OFF."""
import math
import numpy as np
from scipy.optimize import brentq, fsolve

from Avbay_Steady import A, V, NAMES, SIG, Geometry, Settings, board_powers, C2K, CF


def solve_flight(T_amb_c, G_b, G_d, h_o, P):
    """Skin, wall, air and PCB temperatures [K] for one flight condition."""
    geo = Geometry()
    cfg = Settings()
    regions = geo.regions()
    P = np.array(P)
    Aw = geo.inner_area

    T_inf = C2K(T_amb_c)
    T_sky = 0.0552 * T_inf**1.5
    T_rad = T_sky                          # no ground long-wave term in flight
    beta = math.radians(cfg.sun_elevation_deg)

    # absorbed sun: beam on cylinder + diffuse sky only
    q_sun = cfg.solar_absorptivity * (
        G_b * math.cos(beta) / math.pi + G_d / 2.0
    )

    def q_out(Ts):
        return h_o * (Ts - T_inf) + cfg.exterior_emissivity * SIG * (Ts**4 - T_rad**4)

    def skins(Tw):
        return np.array([
            brentq(lambda Ts: ao * (q_out(Ts) - q_sun) - (Tw - Ts) / rw,
                   1.0, max(Tw, T_inf) + 1500.0, xtol=1e-10)
            for _, ao, rw in regions
        ])

    def conduction(Tw, Ts):
        return np.array([(Tw - t) / rw for t, (_, _, rw) in zip(Ts, regions)])

    T_w = brentq(lambda t: conduction(t, skins(t)).sum() - P.sum(),
                 1.0, T_inf + 1500.0, xtol=1e-10)
    T_s = skins(T_w)

    def air(Tp):
        return (cfg.h_board * (A @ Tp) + cfg.h_wall * Aw * T_w) / (
            cfg.h_board * A.sum() + cfg.h_wall * Aw
        )

    def board_loss(Tp):
        w = cfg.board_emissivity * A / Aw
        J_w = SIG * (
            cfg.wall_emissivity * T_w**4
            + (1 - cfg.wall_emissivity) * (w @ Tp**4)
        ) / (
            cfg.wall_emissivity + (1 - cfg.wall_emissivity) * w.sum()
        )
        return (
            cfg.h_board * A * (Tp - air(Tp))
            + cfg.board_emissivity * A * (SIG * Tp**4 - J_w)
        )

    T_p, _, ok, msg = fsolve(lambda Tp: P - board_loss(Tp),
                             np.full(len(P), T_w + 20), xtol=1e-10, full_output=True)
    if ok != 1:
        raise RuntimeError(f'PCB solve failed: {msg}')

    Qconv = np.array([h_o * ao * (t - T_inf) for t, (_, ao, _) in zip(T_s, regions)])
    Qrad = np.array([cfg.exterior_emissivity * SIG * ao * (t**4 - T_rad**4)
                     for t, (_, ao, _) in zip(T_s, regions)])
    Qsolar = np.array([q_sun * ao for _, ao, _ in regions])
    residual = float(P.sum() + Qsolar.sum() - Qconv.sum() - Qrad.sum())
    if abs(residual) > 1e-6:
        raise RuntimeError(f'Energy balance failed: {residual:g} W')

    return dict(P=P, h_o=h_o, T_inf=T_inf, T_rad=T_rad, q_sun=q_sun,
                T_s=T_s[0], T_w=T_w, T_a=air(T_p), T_p=T_p,
                settings=cfg, max_residual_W=abs(residual))


def report(title, r):
    """Flight temperatures and matching ANSYS boundary conditions."""
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
        ('Radiation', 'outer shell', f"e {cfg.exterior_emissivity:g}, sky {CF(r['T_rad'])}"),
        *[('Internal Heat Gen' if i == 0 else '', n, f'{q:,.0f} W/m^3')
          for i, (n, q) in enumerate(zip(NAMES, Q_GEN))],
        ('Convection', 'inside faces', f"h {cfg.h_board:g} W/m^2C, ambient {CF(r['T_a'])}"),
        ('Radiation', 'inside faces', f"Surface to Surface, e {cfg.board_emissivity:g}, wall {CF(r['T_w'])}"),
    ]
    for bc, scope, value in bcs:
        print(f'  {bc:<20}{scope:<14}{value}')


P_PCB = board_powers()

# Flight A: hot / low cooling
T_AMB_C = 40.0
G_B, G_D = 900.0, 100.0
H_O = 25.0
report('FLIGHT A - HOT / LOW COOLING',
       solve_flight(T_AMB_C, G_B, G_D, H_O, P_PCB))

# Flight B: moderate
T_AMB_C = 20.0
G_B, G_D = 1000.0, 100.0
H_O = 75.0
report('FLIGHT B - MODERATE',
       solve_flight(T_AMB_C, G_B, G_D, H_O, P_PCB))

# Flight C: cool / high cooling
T_AMB_C = 0.0
G_B, G_D = 1000.0, 50.0
H_O = 150.0
report('FLIGHT C - COOL / HIGH COOLING',
       solve_flight(T_AMB_C, G_B, G_D, H_O, P_PCB))
