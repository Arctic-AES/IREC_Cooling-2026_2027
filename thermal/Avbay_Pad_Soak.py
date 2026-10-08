"""Pre-launch pad-soak cases for the IREC avionics bay.

Electronics are OFF in these cases. Ground long-wave radiation and
reflected solar are retained, so this file represents environmental heat soak
while the rocket is on the pad.

This script reuses the validated thermal network, geometry, and material
assumptions from Avbay_Steady.py. It intentionally does NOT carry electronics
heat into the pad cases.

Run:
    python thermal/Avbay_Pad_Soak.py
"""

import numpy as np

from Avbay_Steady import CF, NAMES, V, Settings, solve


# name, ambient [C], ground [C], direct normal solar [W/m^2],
# diffuse solar [W/m^2], exterior h [W/(m^2 K)]
PAD_CASES = (
    ("HOT PAD - ELECTRONICS OFF", 42.0, 60.0, 900.0, 100.0, 10.0),
    ("MODERATE PAD - ELECTRONICS OFF", 34.0, 52.0, 900.0, 100.0, 21.6),
    ("COOL / NO SUN - ELECTRONICS OFF", 20.5, 20.5, 0.0, 0.0, 22.0),
)


ELECTRONICS_OFF_W = np.zeros(3)


def report(title, result):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)
    print("Mode: PRE-LAUNCH PAD SOAK")
    print("Electronics: OFF")
    print(f"Ambient:    {CF(result['T_inf'])}")
    print(f"Sky:        {CF(result['T_sky'])}")
    print(f"Radiative surroundings: {CF(result['T_rad'])}")
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
    print("IREC AVIONICS BAY - PAD SOAK")
    print("Ground radiation/reflected solar INCLUDED.")
    print("Electronics heat generation set to ZERO for every pad case.")
    print("These are steady-state environmental soak cases, not transient warm-up times.")

    cfg = Settings()

    for title, T_amb_c, T_gnd_c, G_b, G_d, h_out in PAD_CASES:
        result = solve(
            T_amb_c,
            T_gnd_c,
            G_b,
            G_d,
            h_out,
            P=ELECTRONICS_OFF_W,
            settings=cfg,
        )
        report(title, result)


if __name__ == "__main__":
    main()
