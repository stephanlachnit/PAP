import math as m
import measure as ms

T0 = 273.15
T = 24 + T0
d_T = 1.0

R = 8.3144598
k_Luft = 1.4
M_Luft = 29e-3
c_Luft = m.sqrt(k_Luft * R * T0 / M_Luft)
k_CO2 = 1.3
M_CO2 = 44e-3
c_CO2 = m.sqrt(k_CO2 * R * T0 / M_CO2)

h0_xy = [3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.4e-2]
h0_xy = ms.mean_value(h0_xy)
h0_yt = 3.5e-2
d_h0 = m.sqrt(2/7) * 0.1e-2
c_xy = h0_xy * 10e3
c_yt = h0_yt * 10e3
d_c = d_h0 * 10e3

c0_xy = c_xy * m.sqrt(T0 / T)
d_c0_xy = m.sqrt(T0 / T + (d_c**2 + 0.25 * (c_xy * d_T / T)**2))
sigma_xy = abs(c_Luft - c0_xy) / d_c0_xy
c0_yt = c_yt * m.sqrt(T0 / T)
d_c0_yt = m.sqrt(T0 / T + (d_c**2 + 0.25 * (c_yt * d_T / T)**2))
sigma_yt = abs(c_Luft - c0_yt) / d_c0_yt

print()
print(c_Luft)
print(c_CO2)
print()
print(c_xy)
print(c_yt)
print(d_c)
print()
print(c0_xy)
print(d_c0_xy)
print(sigma_xy)
print()
print(c0_yt)
print(d_c0_yt)
print(sigma_yt)
