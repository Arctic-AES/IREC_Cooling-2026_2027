"""Simple steady ventilation estimate for two 40x10 mm fans."""

from Avbay_Steady import board_powers

# Current sealed hot-case bay-air temperature
T_INITIAL_C = 67.56
T_AMBIENT_C = 42.0

# Two 40x10 mm fans: one intake + one exhaust
FLOW_M3_H = 8.93       # through-flow; fan free-air rating
POWER_PER_FAN_W = 0.30
N_FANS = 2

# Air properties near the hot operating range
RHO_AIR = 1.10         # kg/m^3
CP_AIR = 1007.0        # J/(kg K)

P_ELECTRONICS_W = float(board_powers().sum())
P_FANS_W = N_FANS * POWER_PER_FAN_W
P_TOTAL_W = P_ELECTRONICS_W + P_FANS_W

# Steady ventilation energy balance:
# Q = m_dot * cp * (T_final - T_ambient)
V_DOT = FLOW_M3_H / 3600.0
M_DOT = RHO_AIR * V_DOT
T_FINAL_C = T_AMBIENT_C + P_TOTAL_W / (M_DOT * CP_AIR)
COOLING_C = max(0.0, T_INITIAL_C - T_FINAL_C)

print('2x 40x10 mm fans - 1 intake + 1 exhaust')
print(f'Initial bay air: {T_INITIAL_C:.2f} C')
print(f'Outside ambient: {T_AMBIENT_C:.2f} C')
print(f'Fan heat: {P_FANS_W:.2f} W')
print(f'Flow: {FLOW_M3_H:.2f} m^3/h')
print(f'Predicted ventilated bay air: {T_FINAL_C:.2f} C')
print(f'Cooling from current bay temperature: {COOLING_C:.2f} C')

# Conservative simplifications:
# - electronics remain at full 12.061 W operating heat
# - 100% of fan electrical power becomes heat in the bay
# - rated free-air flow is achieved
# - bay air is perfectly mixed
# - no extra wall/radiation heat rejection is credited
