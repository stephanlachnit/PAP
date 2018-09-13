import math as m
import measure as ms

h0_xy = [3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.5e-2, 3.4e-2]
h0_xy = ms.mean_value(h0_xy)
h0_yt = 3.5e-2
d_h0 = m.sqrt(2/7) * 0.1e-2
c_xy = h0_xy * 10e3
c_yt = h0_yt * 10e3
d_c = d_h0 * 10e3

k = 1.4
R = 8.3144598
T = 24 + 273.15
M = 29e-3
c = m.sqrt(k * R * T / M)

print()
print(c_xy)
print(d_c)
print()
print(c)