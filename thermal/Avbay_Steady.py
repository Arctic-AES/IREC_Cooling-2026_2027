"""
Steady 1D avionics-bay thermal model.
Run: python thermal/Avbay_Steady.py --help
"""
import argparse
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.optimize import brentq, fsolve

SIG = 5.670374419e-8
NAMES = ('C6', 'G6', 'ENIAC')

A = np.array([30145.299, 29958.442, 12339.253]) * 1e-6
V = np.array([21856.975670322, 21875.838207068, 9393.696251272]) * 1e-9

COMPONENT_HEAT_W = {
    'C6': {'U8 TLV1117': 1.659, 'TPS563200 converter-loss allocation': .970,
           'U2 STM32H7A3': .437, 'J4 microSD': .330, 'U1 STM32H562': .286,
           'U6 ADV7280A': .249, 'LP5907-2.8': .137, 'U7 MAX3490': .132,
           'LP5907-1.8': .129, 'Y1 SiT5000': .102, 'Y2 SiT5000': .102},
    'G6': {'R8 24 ohm': .175, 'R9 24 ohm': .175,
           'U1 ISL83488': .131, 'U5 NEO-M9N': .119},
    'ENIAC': {'U5 TL780': 2.619, 'U3 BMT321 (unverified RF-loss estimate)': 1.860,
              'U6 TL780': .622, 'U4 LMX2541': .429, 'U2 HMC361': .415,
              'U7 NCV8163': .290, 'U11 QPL9504': .275},
}
SMALL_PARTS_HEAT_W = {'C6': .194, 'G6': .091, 'ENIAC': .133}


def board_powers():
    return np.array([sum(COMPONENT_HEAT_W[n].values()) + SMALL_PARTS_HEAT_W[n]
                     for n in NAMES])


@dataclass(frozen=True)
class Geometry:
    body_od_m: float = .156718
    body_id_m: float = .152400
    coupler_od_m: float = .152400
    coupler_id_m: float = .147828
    body_length_m: float = .40005
    coupler_length_m: float = .508
    modeled_length_m: float = .40005
    body_k: float = .29
    coupler_k: float = .29
    contact_h: float | None = None

    def regions(self):
        vals = (self.body_od_m, self.body_id_m, self.coupler_od_m,
                self.coupler_id_m, self.body_length_m, self.coupler_length_m,
                self.modeled_length_m, self.body_k, self.coupler_k)
        if not all(math.isfinite(x) and x > 0 for x in vals):
            raise ValueError('Lengths, diameters and conductivities must be positive.')
        if not self.body_od_m > self.body_id_m >= self.coupler_od_m > self.coupler_id_m:
            raise ValueError('Invalid nested tube diameters.')
        if not math.isclose(self.body_id_m, self.coupler_od_m, rel_tol=0, abs_tol=1e-9):
            raise ValueError('A geometric gap needs an explicit gap/adhesive model.')
        if self.modeled_length_m > self.coupler_length_m:
            raise ValueError('Modeled length exceeds the uploaded coupler.')
        if self.contact_h is not None and (not math.isfinite(self.contact_h) or self.contact_h <= 0):
            raise ValueError('Contact conductance must be positive or None.')

        length = min(self.modeled_length_m, self.body_length_m)
        rb = math.log(self.body_od_m / self.body_id_m) / (2 * math.pi * self.body_k * length)
        rc = math.log(self.coupler_od_m / self.coupler_id_m) / (2 * math.pi * self.coupler_k * length)
        ri = 0 if self.contact_h is None else 1 / (self.contact_h * math.pi * self.body_id_m * length)
        regions = [('Body tube + coupler overlap', math.pi * self.body_od_m * length, rb + rc + ri)]

        extension = self.modeled_length_m - length
        if extension > 1e-12:
            rext = math.log(self.coupler_od_m / self.coupler_id_m) / (2 * math.pi * self.coupler_k * extension)
            regions.append(('EXPOSED coupler extension', math.pi * self.coupler_od_m * extension, rext))
        return regions

    @property
    def inner_area(self):
        return math.pi * self.coupler_id_m * self.modeled_length_m


@dataclass(frozen=True)
class Settings:
    sun_elevation_deg: float = 45.0
    ground_albedo: float = .40
    solar_absorptivity: float = .60
    exterior_emissivity: float = .90
    board_emissivity: float = .90
    wall_emissivity: float = .90
    h_board: float = 3.0
    h_wall: float = 3.0
    mount_conductance: tuple = (0., 0., 0.)


def C2K(t): return t + 273.15


def CF(t):
    c = t - 273.15
    return f'{c:.2f} C ({c * 1.8 + 32:.2f} F)'


def solve(T_amb_c, T_gnd_c, G_b, G_d, h_o, P=None, *, geometry=None, settings=None):
    """Skin, wall, air and PCB temperatures [K] for one steady condition."""
    geo, cfg = geometry or Geometry(), settings or Settings()
    regions = geo.regions()
    P = board_powers() if P is None else np.asarray(P, dtype=float)
    mounts = np.asarray(cfg.mount_conductance, dtype=float)

    if P.shape != (3,) or not np.all(np.isfinite(P)) or np.any(P < 0):
        raise ValueError('P must contain three finite nonnegative board loads.')

    Aw = geo.inner_area
    T_inf, T_ground = C2K(T_amb_c), C2K(T_gnd_c)
    T_sky = .0552 * T_inf**1.5
    T_rad = ((T_sky**4 + T_ground**4) / 2) ** .25
    beta = math.radians(cfg.sun_elevation_deg)
    Gh = G_b * math.sin(beta) + G_d
    q_sun = cfg.solar_absorptivity * (
        G_b * math.cos(beta) / math.pi + G_d / 2 + cfg.ground_albedo * Gh / 2
    )

    def q_out(Ts):
        return h_o * (Ts - T_inf) + cfg.exterior_emissivity * SIG * (Ts**4 - T_rad**4)

    def skins(Tw):
        return np.array([
            brentq(lambda Ts: ao * (q_out(Ts) - q_sun) - (Tw - Ts) / rw,
                   1., max(Tw, T_inf, T_ground) + 1500, xtol=1e-10)
            for _, ao, rw in regions
        ])

    def conduction(Tw, Ts):
        return np.array([(Tw - t) / rw for t, (_, _, rw) in zip(Ts, regions)])

    Tw = brentq(lambda t: conduction(t, skins(t)).sum() - P.sum(),
                1., max(T_inf, T_ground) + 1500, xtol=1e-10)
    Ts = skins(Tw)

    def air(Tp):
        return (cfg.h_board * (A @ Tp) + cfg.h_wall * Aw * Tw) / (
            cfg.h_board * A.sum() + cfg.h_wall * Aw
        )

    def board_loss(Tp):
        w = cfg.board_emissivity * A / Aw
        Jw = SIG * (cfg.wall_emissivity * Tw**4 +
                    (1 - cfg.wall_emissivity) * (w @ Tp**4)) / (
                        cfg.wall_emissivity + (1 - cfg.wall_emissivity) * w.sum()
                    )
        conv = cfg.h_board * A * (Tp - air(Tp))
        rad = cfg.board_emissivity * A * (SIG * Tp**4 - Jw)
        mount = mounts * (Tp - Tw)
        return conv + rad + mount

    Tp, _, ok, msg = fsolve(lambda Tp: P - board_loss(Tp),
                            np.full(3, Tw + 20), xtol=1e-10, full_output=True)
    if ok != 1:
        raise RuntimeError(f'PCB solve failed: {msg}')

    Ta = air(Tp)
    Qwall = conduction(Tw, Ts)
    Qconv = np.array([h_o * ao * (t - T_inf) for t, (_, ao, _) in zip(Ts, regions)])
    Qrad = np.array([cfg.exterior_emissivity * SIG * ao * (t**4 - T_rad**4)
                     for t, (_, ao, _) in zip(Ts, regions)])
    Qsolar = np.array([q_sun * ao for _, ao, _ in regions])

    # residual check
    residual = float(P.sum() + Qsolar.sum() - Qconv.sum() - Qrad.sum())
    if abs(residual) > 1e-6:
        raise RuntimeError(f'Energy balance failed: {residual:g} W')

    return dict(P=P, h_o=h_o, T_inf=T_inf, T_rad=T_rad, q_sun=q_sun,
                T_s=Ts[0], T_w=Tw, T_a=Ta, T_p=Tp,
                geometry=geo, settings=cfg, max_residual_W=abs(residual))


def report(title, r):
    """Temperatures and matching ANSYS steady-state boundary conditions."""
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


# Current full-power steady cases
P_PCB = board_powers()

# Hot pad
T_AMB_C = 42.0
T_GND_C = 60.0
G_B, G_D = 900.0, 100.0
H_O = 10.0
HOT_CASE = ('HOT PAD', T_AMB_C, T_GND_C, G_B, G_D, H_O)

# Moderate pad
T_AMB_C = 34.0
T_GND_C = 52.0
G_B, G_D = 900.0, 100.0
H_O = 21.6
MODERATE_CASE = ('MODERATE PAD', T_AMB_C, T_GND_C, G_B, G_D, H_O)

# Cool / no sun
T_AMB_C = 20.5
T_GND_C = 20.5
G_B, G_D = 0.0, 0.0
H_O = 22.0
COOL_CASE = ('COOL / NO SUN', T_AMB_C, T_GND_C, G_B, G_D, H_O)

CASES = (HOT_CASE, MODERATE_CASE, COOL_CASE)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--length-mm', type=float, default=400.05)
    ap.add_argument('--contact-h', type=float, default=None)
    ap.add_argument('--h-inside', type=float, default=3.)
    args = ap.parse_args()

    geo = replace(Geometry(), modeled_length_m=args.length_mm / 1000,
                  contact_h=args.contact_h)
    cfg = replace(Settings(), h_board=args.h_inside, h_wall=args.h_inside)

    for title, T_amb_c, T_gnd_c, G_b, G_d, h_o in CASES:
        report(title, solve(T_amb_c, T_gnd_c, G_b, G_d, h_o, P_PCB,
                            geometry=geo, settings=cfg))


if __name__ == '__main__':
    main()
