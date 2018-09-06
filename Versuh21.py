import math as m
#import measure as ms

# CuSO4
I = 1.0
d_I = 0.03
t = 1898.0
d_t = 3.0
mv = [93.7898e-3, 92.4716e-3]
mn = [93.1120e-3, 93.0984e-3]
d_m = 0.0010e-3
m0 = [mv[0] - mn[0], mn[1] - mv[1]]
d_m0 = m.sqrt(2) * d_m
M = 63.546e-3
d_M = 0.003e-3
z = 2.0
F = [0.0, 0.0]
d_F = [0.0, 0.0]
for i in range(2):
 F[i] = I * t * M / (z * m0[i])
 d_F[i] = m.sqrt(M**2 * ((t * d_I)**2 + (I * d_t)**2 + (I * t * d_m0 / m0[i])**2) + (I * t * d_M)**2) / (z * m0[i])

print()
print(m0)
print(d_m0)
print()
print(F)
print(d_F)

# H2SO4
I = 0.72
d_I = 0.02
t = 438.0
d_t = 3.0
TR = 273.15 + 24.5
d_TR = 1.0
TW = 273.15 + 26.5
d_TW = 1.0
T = 0.5 * (TR + TW)
d_T = T - TR
Vv = [2.0e-6, 3.2e-6]
Vn = [23.2e-6, 46.6e-6]
V = [Vn[0] - Vv[0] , Vn[1] - Vv[1]]
d_V = m.sqrt(2) * 0.3e-6
pL = 1002.5e2
d_pL = 1.0e2
pD = 24.5 * 133.322
d_pD = 1.5 * 133.322
p0 = 1013.25e2
V0 = 22.414e-3
T0 = 273.15
z = [4.0, 2.0]
F = [0.0, 0.0]
d_F = [0.0, 0.0]
for i in range(2):
 F[i] = I * t * p0 * T * V0 / (z[i] * (pL - 0.9 * pD) * T0 * V[i])
 d_F[i] = p0 * V0 / (z[i] * T0) * m.sqrt(1.0 / (pL - 0.9 * pD)**2 * ((t * T / V[i] * d_I)**2 + (I * T / V[i] * d_t)**2 + (I * t / V[i] * d_T)**2 + (I * t * T / V[i]**2 * d_V)**2) + (I * t * T)**2 / (V[i]**2 * (pL - 0.9 * pD)**4) * (d_pL**2 + 0.9**2 * d_pD**2))

print()
print(T)
print(d_T)
print()
print(V)
print(d_V)
print()
print(pD)
print(d_pD)
print()
print(F)
print(d_F)
