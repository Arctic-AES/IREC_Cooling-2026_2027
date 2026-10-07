"""Sealed-bay steady thermal network: ONE isothermal node and heat load per PCB.

Component estimates are retained for bookkeeping and summed once per board.
Predictions are average board temperatures, not chip junction temperatures.
Run: python thermal/Avbay_Steady.py --help. Python 3.10+, numpy, scipy.
"""
import argparse
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.optimize import brentq, fsolve

SIG = 5.670374419e-8
NAMES = ('C6', 'G6', 'ENIAC')
# Bare-board STEP surface areas, including edges; all block bodies omitted.
A = np.array([30145.299, 29958.442, 12339.253]) * 1e-6
V = np.array([21856.975670322, 21875.838207068, 9393.696251272]) * 1e-9
# Preliminary user-provided operating-case heat estimates [W], NOT measured.
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
    """Single source of truth; these totals are the only electronics loads."""
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
    # Provisional analysis domain: common overlap, NOT confirmed sealed length.
    modeled_length_m: float = .40005
    body_k: float = .29
    coupler_k: float = .29
    contact_h: float | None = None  # W/(m^2 K); None means ideal contact

    def regions(self):
        """Return (name, exterior area [m^2], radial resistance [K/W])."""
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
    sun_elevation_deg: float = 45.0  # vertical rocket, uniform circumferential average
    ground_albedo: float = .40
    solar_absorptivity: float = .60
    exterior_emissivity: float = .90
    board_emissivity: float = .90
    wall_emissivity: float = .90
    h_board: float = 3.0
    h_wall: float = 3.0
    mount_conductance: tuple = (0., 0., 0.)  # PCB-to-wall W/K; mounts absent from STEP


def C2K(t): return t + 273.15


def CF(t):
    c = t - 273.15
    return f'{c:.2f} C ({c * 1.8 + 32:.2f} F)'


def solve(T_amb_c, T_gnd_c, G_b, G_d, h_o, P=None, *, geometry=None, settings=None):
    """Passive sealed equilibrium. G_b=direct NORMAL, G_d=diffuse horizontal.

    All returned temperatures are kelvin. One inner-wall and one air node.
    Optional extension is assumed exposed to the same environment as the tube.
    """
    geo, cfg = geometry or Geometry(), settings or Settings()
    regions = geo.regions()
    P = board_powers() if P is None else np.asarray(P, dtype=float)
    mounts = np.asarray(cfg.mount_conductance, dtype=float)
    if P.shape != (3,) or not np.all(np.isfinite(P)) or np.any(P < 0):
        raise ValueError('P must contain three finite nonnegative board loads.')
    if mounts.shape != (3,) or not np.all(np.isfinite(mounts)) or np.any(mounts < 0):
        raise ValueError('Mount conductances must be three finite nonnegative values.')
    if not all(math.isfinite(x) and x >= 0 for x in (G_b, G_d, h_o)):
        raise ValueError('Solar irradiances and exterior h must be finite/nonnegative.')
    if not all(math.isfinite(x) for x in (T_amb_c, T_gnd_c)):
        raise ValueError('Environmental temperatures must be finite.')
    for e in (cfg.exterior_emissivity, cfg.board_emissivity, cfg.wall_emissivity):
        if not math.isfinite(e) or not 0 < e <= 1:
            raise ValueError('Emissivities must lie in (0, 1].')
    for f in (cfg.ground_albedo, cfg.solar_absorptivity):
        if not math.isfinite(f) or not 0 <= f <= 1:
            raise ValueError('Albedo and solar absorptivity must lie in [0, 1].')
    if not math.isfinite(cfg.sun_elevation_deg) or not 0 <= cfg.sun_elevation_deg <= 90:
        raise ValueError('Sun elevation must lie in [0, 90] degrees.')
    if not all(math.isfinite(h) and h > 0 for h in (cfg.h_board, cfg.h_wall)):
        raise ValueError('Internal convection coefficients must be positive.')
    Aw = geo.inner_area
    if A.sum() > Aw:
        raise ValueError('Board-to-wall view-factor approximation invalid: board area exceeds wall area.')
    T_inf, T_ground = C2K(T_amb_c), C2K(T_gnd_c)
    if min(T_inf, T_ground) <= 0:
        raise ValueError('Temperatures must exceed absolute zero.')
    T_sky = .0552 * T_inf**1.5  # clear-sky Swinbank approximation; kelvin
    T_rad = ((T_sky**4 + T_ground**4) / 2) ** .25
    beta = math.radians(cfg.sun_elevation_deg)
    Gh = G_b * math.sin(beta) + G_d
    q_sun = cfg.solar_absorptivity * (G_b * math.cos(beta) / math.pi + G_d / 2 + cfg.ground_albedo * Gh / 2)

    def q_out(Ts):
        return h_o * (Ts - T_inf) + cfg.exterior_emissivity * SIG * (Ts**4 - T_rad**4)

    def skins(Tw):
        return np.array([brentq(lambda Ts: ao * (q_out(Ts) - q_sun) - (Tw - Ts) / rw,
                                1., max(Tw, T_inf, T_ground) + 1500, xtol=1e-10)
                         for _, ao, rw in regions])

    def conduction(Tw, Ts):
        return np.array([(Tw - t) / rw for t, (_, _, rw) in zip(Ts, regions)])

    # Total net sidewall conduction = electronics heat. Solar enters exterior nodes.
    Tw = brentq(lambda t: conduction(t, skins(t)).sum() - P.sum(),
                1., max(T_inf, T_ground) + 1500, xtol=1e-10)
    Ts = skins(Tw)

    def air(Tp):
        return (cfg.h_board * (A @ Tp) + cfg.h_wall * Aw * Tw) / (cfg.h_board * A.sum() + cfg.h_wall * Aw)

    def terms(Tp):
        # Board-to-wall Fi,w=1; reciprocal Fw,i=Ai/Aw. Board-to-board radiation omitted.
        w = cfg.board_emissivity * A / Aw
        Jw = SIG * (cfg.wall_emissivity * Tw**4 + (1 - cfg.wall_emissivity) * (w @ Tp**4)) / (
            cfg.wall_emissivity + (1 - cfg.wall_emissivity) * w.sum())
        conv = cfg.h_board * A * (Tp - air(Tp))
        rad = cfg.board_emissivity * A * (SIG * Tp**4 - Jw)
        mount = mounts * (Tp - Tw)
        return conv, rad, mount

    Tp, _, flag, msg = fsolve(lambda t: sum(terms(t)) - P, np.full(3, Tw + 20), xtol=1e-10, full_output=True)
    conv, rad, mount = terms(Tp)
    Ta = air(Tp)
    Qwall = conduction(Tw, Ts)
    Qconv = np.array([h_o * ao * (t - T_inf) for t, (_, ao, _) in zip(Ts, regions)])
    Qrad = np.array([cfg.exterior_emissivity * SIG * ao * (t**4 - T_rad**4) for t, (_, ao, _) in zip(Ts, regions)])
    Qsolar = np.array([q_sun * ao for _, ao, _ in regions])
    residuals = {
        'boards_W': P - conv - rad - mount,
        'air_W': float(conv.sum() - cfg.h_wall * Aw * (Ta - Tw)),
        'inner_wall_W': float(cfg.h_wall * Aw * (Ta - Tw) + rad.sum() + mount.sum() - Qwall.sum()),
        'outer_regions_W': Qwall + Qsolar - Qconv - Qrad,
        'whole_bay_W': float(P.sum() + Qsolar.sum() - Qconv.sum() - Qrad.sum()),
    }
    err = max(float(np.max(np.abs(x))) for x in residuals.values())
    if not np.all(np.isfinite(Tp)) or np.any(Tp <= 0) or not math.isfinite(err) or err > 1e-6:
        raise RuntimeError(f'Thermal solve failed (flag {flag}): {msg}; residual {err:g} W')
    return dict(P=P, h_o=h_o, T_inf=T_inf, T_rad=T_rad, T_sky=T_sky, q_sun=q_sun,
                T_s=Ts[0], T_s_regions=Ts, T_w=Tw, T_a=Ta, T_p=Tp,
                geometry=geo, settings=cfg, regions=regions, residuals=residuals,
                max_residual_W=err, Q_solar_W=Qsolar.sum(),
                Q_ext_convection_W=Qconv.sum(), Q_ext_radiation_W=Qrad.sum())


# Same electronics operating-case budget in every environmental scenario.
# Idle-mode heat must be audited separately, not inferred from voltage scaling.
CASES = (
    ('HOT PAD (assumed scenario, not a verified worst-case bound)', 42., 60., 900., 100., 10.),
    ('MODERATE PAD', 34., 52., 900., 100., 21.6),
    ('COOL / NO SUN', 20.5, 20.5, 0., 0., 22.),
)


def report(title, r):
    geo, cfg = r['geometry'], r['settings']
    print(f"\n{title}\nAmbient: {CF(r['T_inf'])}\nBay air: {CF(r['T_a'])}\nInner wall: {CF(r['T_w'])}")
    for (label, ao, rw), t in zip(r['regions'], r['T_s_regions']):
        print(f'{label}: skin {CF(t)}, outside area {ao:.6f} m^2, radial R {rw:.6f} K/W')
    print(f'Inner cylindrical area: {geo.inner_area:.6f} m^2')
    print('PCB      Total heat [W]  FULL-board IHG [W/m^3]  Average PCB temperature')
    for n, p, v, t in zip(NAMES, r['P'], V, r['T_p']):
        print(f'{n:<8}{p:14.3f}{p/v:24,.0f}  {CF(t)}')
    print(f"Electronics {r['P'].sum():.3f} W; solar {r['Q_solar_W']:.3f} W")
    print(f"Exterior convection {r['Q_ext_convection_W']:.3f} W; radiation {r['Q_ext_radiation_W']:.3f} W (positive outward)")
    print(f"Maximum energy-balance residual: {r['max_residual_W']:.3e} W")
    print('ANSYS BASELINE: simplified boards, no active component loads')
    print(f"  Exterior cylindrical faces: inward solar Heat Flux {r['q_sun']:.3f} W/m^2")
    print(f"  Outside Convection: h={r['h_o']:g} W/(m^2 K), bulk={CF(r['T_inf'])}")
    print(f"  Outside Radiation / To Ambient: e={cfg.exterior_emissivity:g}, surroundings={CF(r['T_rad'])}")
    print('  FULL-board IHG from table, once per PCB; suppress component IHG and Heat Flow.')
    print(f"  Internal Convection: board h={cfg.h_board:g}, wall h={cfg.h_wall:g}; start bulk at {CF(r['T_a'])}.")
    print('  Iterate bulk=sum(h*A*Tmean)/sum(h*A) until total internal convection is zero.')
    print(f'  Internal surface radiation: PCB e={cfg.board_emissivity:g}, wall e={cfg.wall_emissivity:g}.')
    print('  Use actual view factors; Perfect only after including real closed bay boundaries.')
    print('  Do not prescribe inner-wall/board temperatures; initial temperature does not set steady results.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--length-mm', type=float, default=400.05,
                    help='Provisional modeled length; above 400.05 assumes EXPOSED coupler extension.')
    ap.add_argument('--contact-h', type=float, default=None, help='Tube/coupler contact W/(m^2 K); default ideal')
    ap.add_argument('--h-inside', type=float, default=3., help='Assumed internal h W/(m^2 K)')
    args = ap.parse_args()
    geo = replace(Geometry(), modeled_length_m=args.length_mm / 1000, contact_h=args.contact_h)
    cfg = replace(Settings(), h_board=args.h_inside, h_wall=args.h_inside)
    geo.regions()
    print('PRELIMINARY: average board temperatures, no component hot spots/junction predictions.')
    print(f'Domain length {args.length_mm:g} mm; real sealed length and bulkheads not confirmed.')
    print('Endcap heat transfer and mounts omitted. Tube contact ideal unless --contact-h supplied.')
    print('Coupler extension is EXPOSED.' if args.length_mm > 400.05 else '108 mm STEP coupler extension is outside default domain.')
    print('COMPONENT BUDGET: estimates only, NOT extra simulation loads')
    for n in NAMES:
        print(f'\n{n}')
        for label, p in COMPONENT_HEAT_W[n].items(): print(f'  {label:<54}{p:.3f} W')
        print(f'  Remaining small parts: {SMALL_PARTS_HEAT_W[n]:.3f} W')
    print('\nSame component budget in all ambient scenarios; idle-mode power not yet verified.')
    for title, *env in CASES:
        report(title, solve(*env, geometry=geo, settings=cfg))


if __name__ == '__main__':
    main()
