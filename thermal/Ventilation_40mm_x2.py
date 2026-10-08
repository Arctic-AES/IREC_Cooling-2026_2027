"""Vent-hole sensitivity for two 40x10 mm 5 V fans."""

import math
from scipy.optimize import brentq
from Avbay_Steady import board_powers

T_INITIAL_C = 67.56
T_AMB_C = 42.0

RHO_AIR = 1.10
CP_AIR = 1007.0
CD = 0.62

FLOW_FREE_M3_H = 8.93
STATIC_PER_FAN_MMH2O = 1.95
POWER_PER_FAN_W = 0.30
N_FANS = 2

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
    q_free = FLOW_FREE_M3_H / 3600.0
    p0_pair = N_FANS * STATIC_PER_FAN_MMH2O * 9.80665
    return p0_pair * max(0.0, 1.0 - q_m3_s / q_free)


def vent_pressure_pa(q_m3_s, area_in, area_out):
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
    p_total = float(board_powers().sum()) + N_FANS * POWER_PER_FAN_W
    mdot = RHO_AIR * (flow_m3_h / 3600.0)
    return T_AMB_C + p_total / (mdot * CP_AIR)


print('2x 40x10 mm fans: opposite-side inlet + outlet')
print(f'Initial bay air: {T_INITIAL_C:.2f} C | Ambient: {T_AMB_C:.2f} C')
print(f'Fan heat: {N_FANS * POWER_PER_FAN_W:.2f} W')
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
