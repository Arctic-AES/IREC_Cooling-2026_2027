"""Pre-launch pad-soak cases for the IREC avionics bay.

Electronics are OFF in these cases. Ground long-wave radiation and
reflected solar are retained, so this file represents environmental heat soak
while the rocket is on the pad.

Run:
    python thermal/Avbay_Pad_Soak.py
"""

import numpy as np

from Avbay_Steady import Settings, solve


PAD_CASES = (
    ("HOT PAD - ELECTRONICS OFF", 42.0, 60.0, 900.0, 100.0, 10.0),
    ("MODERATE PAD - ELECTRONICS OFF", 34.0, 52.0, 900.0, 100.0, 21.6),
    ("COOL / NO SUN - ELECTRONICS OFF", 20.5, 20.5, 0.0, 0.0, 22.0),
)

ELECTRONICS_OFF_W = np.zeros(3)


def report(title, result, T_gnd_c, h_out):
    print(f"\n{title}")
    print(
        f"Tamb={result['T_inf'] - 273.15:.2f} C | "
        f"Tg={T_gnd_c:.2f} C | "
        f"h={h_out:.2f} W/m2K | "
        f"qsolar={result['q_sun']:.2f} W/m2"
    )
    print(
        f"Ts={result['T_s'] - 273.15:.2f} C | "
        f"Tw={result['T_w'] - 273.15:.2f} C | "
        f"Ta={result['T_a'] - 273.15:.2f} C"
    )
    print(
        f"C6={result['T_p'][0] - 273.15:.2f} C | "
        f"G6={result['T_p'][1] - 273.15:.2f} C | "
        f"ENIAC={result['T_p'][2] - 273.15:.2f} C"
    )
    print(
        f"P={result['P'].sum():.3f} W | "
        f"residual={result['max_residual_W']:.2e} W"
    )


def main():
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
        report(title, result, T_gnd_c, h_out)


if __name__ == "__main__":
    main()
