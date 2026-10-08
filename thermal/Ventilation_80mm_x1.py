"""Vent-hole sensitivity for one 80x25 mm 5 V fan.

The fan drives outside air through an inlet bank on one side of the bay and
an outlet bank on the opposite side. This is side-to-side ventilation, not
top-to-bottom natural draft.
"""

import math
from scipy.optimize import brentq
from Avbay_Steady import board_powers

T_INITIAL_C = 67.56
T_AMB_C = 42.0

RHO_AIR = 1.10      # kg/m^3
CP_AIR = 1007.0     # J/(kg K)
CD = 0.62           # sharp-edged round-hole assumption

# Noctua NF-A8 5V class data
FLOW_FREE_M3_H = 55.5
STATIC_MMH2O = 2.37
FAN_POWER_W = 0.74

# Identical inlet and outlet banks on opposite sides of the bay.
# (holes per side, hole diameter mm)
VENT_STAGES = (
    (8, 6),
    (8, 8),
    (12, 8),
    (12, 10),
    (16, 10),
    (20, 10),
    (20, 12),
)


def open_area(n_holes, diameter_mm):
    d = diameter_mm / 1000.0
    return n_holes * math.pi * d**2 / 4.0


def fan_pressure_pa(q_m3_s):
    """Linearized fan curve from shutoff pressure to free-air flow."""
    q_free = FLOW_FREE_M3_H / 3600.0
    p0 = STATIC_MMH2O * 9.80665
    return p0 * max(0.0, 1.0 - q_m3_s / q_free)


def vent_pressure_pa(q_m3_s, area_in, area_out):
    """Sharp-edged inlet + outlet losses in series."""
    vin = q_m3_s / (CD * area_in)
    vout = q_m3_s / (CD * area_out)
    return 0.5 * RHO_AIR * (vin**2 + vout**2)


def actual_flow_m3_h(n_holes, diameter_mm):
    area = open_area(n_holes, diameter_mm)
    q_free = FLOW_FREE_M3_H / 3600.0
    q = brentq(
        lambda x: fan_pressure_pa(x) - vent_pressure_pa(x, area, area),
        0.0,
        q_free,
    )
    return q * 3600.0


def final_air_temp(flow_m3_h):
    """Q = m_dot cp (Tbay - Tamb), including fan electrical heat."""
    p_total = float(board_powers().sum()) + FAN_POWER_W
    mdot = RHO_AIR * (flow_m3_h / 3600.0)
    return T_AMB_C + p_total / (mdot * CP_AIR)


print('1x 80x25 mm fan: opposite-side inlet + outlet')
print(f'Initial bay air: {T_INITIAL_C:.2f} C | Ambient: {T_AMB_C:.2f} C')
print(f'Fan heat: {FAN_POWER_W:.2f} W')
print('holes/side  dia   area/side     flow       final air   cooling')

for n, d in VENT_STAGES:
    area_mm2 = open_area(n, d) * 1e6
    flow = actual_flow_m3_h(n, d)
    final_t = final_air_temp(flow)
    cooling = T_INITIAL_C - final_t
    print(
        f'{n:>5}      {d:>2} mm   {area_mm2:>7.0f} mm^2   '
        f'{flow:>6.2f} m^3/h   {final_t:>6.2f} C   {cooling:>6.2f} C'
    )

print('\nDTEG 8.1 requires adequate venting and VFRR discussion; it does not set a specific hole diameter/count.')
print('Flow is a first-pass estimate using a linear fan curve and sharp-edged-hole losses; test the final vent geometry.')
