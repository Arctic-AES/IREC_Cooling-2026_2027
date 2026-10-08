"""In-flight steady thermal sensitivity cases for the IREC avionics bay.

Electronics are ON in these cases. Pad-soak ground heating is REMOVED:
  * no ground long-wave radiation term
  * no ground-reflected solar term

The flight cases are sensitivity points, not a full trajectory model.
Exterior convection coefficient and ambient temperature should later be
replaced by values derived from RocketPy/trajectory data and/or Fluent.

Run:
    python thermal/Avbay_Flight.py
"""

import math
import numpy as np
from scipy.optimize import brentq, fsolve

from Avbay_Steady import (
    A,
    V,
    NAMES,
    SIG,
    Geometry,
    Settings,
    board_powers,
    C2K,
    CF,
)


# name, ambient [C], direct normal solar [W/m^2], diffuse solar [W/m^2],
# exterior h [W/(m^2 K)]
# These are sensitivity points until replaced by actual trajectory-derived data.
FLIGHT_CASES = (
    ("FLIGHT A - HOT / LOW COOLING", 40.0, 900.0, 100.0, 25.0),
    ("FLIGHT B - MODERATE", 20.0, 1000.0, 100.0, 75.0),
    ("FLIGHT C - COOL / HIGH COOLING", 0.0, 1000.0, 50.0, 150.0),
)


def solve_flight(T_amb_c, G_b, G_d, h_out, *, geometry=None, settings=None):
    """Steady in-flight thermal solve with electronics ON and no ground term."""
    geo = geometry or Geometry()
    cfg = settings or Settings()
    regions = geo.regions()
    P = board_powers()
    mounts = np.asarray(cfg.mount_conductance, dtype=float)
    Aw = geo.inner_area

    T_inf = C2K(T_amb_c)
    T_sky = 0.0552 * T_inf**1.5
    T_rad = T_sky

    beta = math.radians(cfg.sun_elevation_deg)

    # In flight: direct beam averaged around vertical cylinder + diffuse sky only.
    # No ground-reflected solar term.
    q_sun = cfg.solar_absorptivity * (
        G_b * math.cos(beta) / math.pi
        + G_d / 2.0
    )

    def q_out(Ts):
        return (
            h_out * (Ts - T_inf)
            + cfg.exterior_emissivity * SIG * (Ts**4 - T_rad**4)
        )

    def skins(Tw):
        upper = max(Tw, T_inf) + 1500.0
        return np.array([
            brentq(
                lambda Ts: ao * (q_out(Ts) - q_sun) - (Tw - Ts) / rw,
                1.0,
                upper,
                xtol=1e-10,
            )
            for _, ao, rw in regions
        ])

    def conduction(Tw, Ts):
        return np.array([(Tw - t) / rw for t, (_, _, rw) in zip(Ts, regions)])

    Tw = brentq(
        lambda t: conduction(t, skins(t)).sum() - P.sum(),
        1.0,
        T_inf + 1500.0,
        xtol=1e-10,
    )
    Ts = skins(Tw)

    def air(Tp):
        return (
            cfg.h_board * (A @ Tp) + cfg.h_wall * Aw * Tw
        ) / (
            cfg.h_board * A.sum() + cfg.h_wall * Aw
        )

    def terms(Tp):
        w = cfg.board_emissivity * A / Aw
        Jw = SIG * (
            cfg.wall_emissivity * Tw**4
            + (1.0 - cfg.wall_emissivity) * (w @ Tp**4)
        ) / (
            cfg.wall_emissivity
            + (1.0 - cfg.wall_emissivity) * w.sum()
        )
        conv = cfg.h_board * A * (Tp - air(Tp))
        rad = cfg.board_emissivity * A * (SIG * Tp**4 - Jw)
        mount = mounts * (Tp - Tw)
        return conv, rad, mount

    Tp, _, flag, msg = fsolve(
        lambda t: sum(terms(t)) - P,
        np.full(3, Tw + 20.0),
        xtol=1e-10,
        full_output=True,
    )

    conv, rad, mount = terms(Tp)
    Ta = air(Tp)
    Qwall = conduction(Tw, Ts)
    Qconv = np.array([
        h_out * ao * (t - T_inf)
        for t, (_, ao, _) in zip(Ts, regions)
    ])
    Qrad = np.array([
        cfg.exterior_emissivity * SIG * ao * (t**4 - T_rad**4)
        for t, (_, ao, _) in zip(Ts, regions)
    ])
    Qsolar = np.array([q_sun * ao for _, ao, _ in regions])

    residuals = {
        "boards_W": P - conv - rad - mount,
        "air_W": float(conv.sum() - cfg.h_wall * Aw * (Ta - Tw)),
        "inner_wall_W": float(
            cfg.h_wall * Aw * (Ta - Tw)
            + rad.sum()
            + mount.sum()
            - Qwall.sum()
        ),
        "outer_regions_W": Qwall + Qsolar - Qconv - Qrad,
        "whole_bay_W": float(P.sum() + Qsolar.sum() - Qconv.sum() - Qrad.sum()),
    }
    err = max(float(np.max(np.abs(x))) for x in residuals.values())

    if (
        flag != 1
        or not np.all(np.isfinite(Tp))
        or np.any(Tp <= 0)
        or err > 1e-6
    ):
        raise RuntimeError(
            f"Flight thermal solve failed (flag {flag}): {msg}; residual {err:g} W"
        )

    return {
        "P": P,
        "T_inf": T_inf,
        "T_sky": T_sky,
        "T_rad": T_rad,
        "q_sun": q_sun,
        "T_s_regions": Ts,
        "T_w": Tw,
        "T_a": Ta,
        "T_p": Tp,
        "Q_solar_W": float(Qsolar.sum()),
        "Q_ext_convection_W": float(Qconv.sum()),
        "Q_ext_radiation_W": float(Qrad.sum()),
        "max_residual_W": err,
        "regions": regions,
    }


def report(title, result, h_out):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)
    print("Mode: IN FLIGHT")
    print("Electronics: ON")
    print("Ground radiation/reflection: REMOVED")
    print(f"Ambient:    {CF(result['T_inf'])}")
    print(f"Sky:        {CF(result['T_sky'])}")
    print(f"Exterior h: {h_out:.2f} W/(m^2 K)")
    print(f"Absorbed solar flux: {result['q_sun']:.2f} W/m^2")
    print(f"Bay air:    {CF(result['T_a'])}")
    print(f"Inner wall: {CF(result['T_w'])}")

    for (label, area, resistance), temp in zip(result['regions'], result['T_s_regions']):
        print(
            f"{label}: {CF(temp)} | "
            f"outside area={area:.6f} m^2 | R={resistance:.6f} K/W"
        )

    print("\nPCB      Heat [W]   IHG [W/m^3]   Average PCB temperature")
    for name, power, volume, temp in zip(NAMES, result['P'], V, result['T_p']):
        print(f"{name:<8}{power:9.3f}{power / volume:15,.0f}   {CF(temp)}")

    print(f"Total electronics heat: {result['P'].sum():.3f} W")
    print(f"Solar absorbed:          {result['Q_solar_W']:.3f} W")
    print(f"Exterior convection:     {result['Q_ext_convection_W']:.3f} W")
    print(f"Exterior radiation:      {result['Q_ext_radiation_W']:.3f} W")
    print(f"Max energy residual:     {result['max_residual_W']:.3e} W")


def main():
    print("IREC AVIONICS BAY - FLIGHT THERMAL SENSITIVITY")
    print("Electronics heat generation ON using Avbay_Steady.py board power budget.")
    print("Ground long-wave radiation and reflected solar are removed.")
    print("These are steady sensitivity points, not a transient ascent trajectory.")
    print("Replace ambient/h_out with trajectory or Fluent-derived values when available.")

    for title, T_amb_c, G_b, G_d, h_out in FLIGHT_CASES:
        result = solve_flight(T_amb_c, G_b, G_d, h_out)
        report(title, result, h_out)


if __name__ == "__main__":
    main()
