"""Simple ambient-air ventilation estimate for one 80x25 mm fan."""

import math
from Avbay_Steady import Geometry, board_powers

# Fan data supplied
FLOW_M3_H = 55.5
FAN_POWER_W = 0.74

# Simple hot-pad comparison point
T_AMB_C = 42.0
RHO_AIR = 1.10       # kg/m^3, approximate hot-air density
CP_AIR = 1007.0      # J/(kg K)


def temperature_rise(flow_m3_h, heat_w):
    """Well-mixed steady ventilation: Q = m_dot * cp * delta_T."""
    vdot = flow_m3_h / 3600.0
    mdot = RHO_AIR * vdot
    return heat_w / (mdot * CP_AIR)


def air_change_time(flow_m3_h):
    """Ideal time to move one empty-bay volume at rated free-air flow."""
    geo = Geometry()
    bay_volume = math.pi * (geo.coupler_id_m / 2.0) ** 2 * geo.modeled_length_m
    return bay_volume / (flow_m3_h / 3600.0)


P_ELECTRONICS = float(board_powers().sum())
P_TOTAL = P_ELECTRONICS + FAN_POWER_W
DT = temperature_rise(FLOW_M3_H, P_TOTAL)

print('1x 80x25 mm fan')
print(f'Electronics heat: {P_ELECTRONICS:.3f} W | fan heat: {FAN_POWER_W:.3f} W')
print(f'Ideal flow: {FLOW_M3_H:.2f} m^3/h | dT={DT:.2f} C | air={T_AMB_C + DT:.2f} C | one-volume={air_change_time(FLOW_M3_H):.2f} s')
print('Rated free-air flow is an upper bound; vents and restrictions reduce real flow.')
