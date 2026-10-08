"""Simple ambient-air ventilation estimate for two 40x10 mm fans."""

import math
from Avbay_Steady import Geometry, board_powers

# Fan data supplied for one fan
FLOW_PER_FAN_M3_H = 8.93
POWER_PER_FAN_W = 0.30
N_FANS = 2

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
P_FANS = N_FANS * POWER_PER_FAN_W
P_TOTAL = P_ELECTRONICS + P_FANS

# If both fans move air in parallel with adequate passive inlet/outlet area,
# their free-air ratings can approximately add.
FLOW_PARALLEL = N_FANS * FLOW_PER_FAN_M3_H
DT_PARALLEL = temperature_rise(FLOW_PARALLEL, P_TOTAL)

# If one fan is intake and one is exhaust, they are effectively in series
# through the bay. Their free-air flow ratings do NOT add. Without fan curves,
# use one-fan free-air flow as a simple nominal comparison.
FLOW_IN_OUT = FLOW_PER_FAN_M3_H
DT_IN_OUT = temperature_rise(FLOW_IN_OUT, P_TOTAL)

print('2x 40x10 mm fans')
print(f'Electronics heat: {P_ELECTRONICS:.3f} W | fan heat: {P_FANS:.3f} W')
print(f'Parallel-flow ideal: {FLOW_PARALLEL:.2f} m^3/h | dT={DT_PARALLEL:.2f} C | air={T_AMB_C + DT_PARALLEL:.2f} C | one-volume={air_change_time(FLOW_PARALLEL):.2f} s')
print(f'1 intake + 1 exhaust: ~{FLOW_IN_OUT:.2f} m^3/h | dT={DT_IN_OUT:.2f} C | air={T_AMB_C + DT_IN_OUT:.2f} C | one-volume={air_change_time(FLOW_IN_OUT):.2f} s')
print('Rated free-air flow is an upper bound; vents and restrictions reduce real flow.')
