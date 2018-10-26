import math as m
import measure as ms

eta0 = 1.81e-5
rhoO = 8.704e2
rhoL = 1.29
rho = rhoO - rhoL
g = 9.81
d = 6.00e-3
U = 501
b = 7.78e-3
p = 9.98e4

sf = [10, 10, 10, 10, 10]
tf = [11.64, 13.26, 12.90, 12.33, 12.19]
ss = [10, 10, 10, 10, 10]
ts = [13.92, 16.16, 15.58, 16.59, 15.88]

r = [0.0, 0.0, 0.0, 0.0, 0.0]
q = [0.0, 0.0, 0.0, 0.0, 0.0]
eta = [0.0, 0.0, 0.0, 0.0, 0.0]

for i in range(5):
 sf[i] *= 5e-5
 ss[i] *= 5e-5
 r[i] = 3 * m.sqrt((eta0 * sf[i]) / (2 * rho * g * tf[i]))
 eta[i] = eta0 / (1.0 + b / (r[i] * p))
 q[i] = (sf[i] / tf[i] + ss[i] / ts[i]) * 18 * m.pi * d / U * m.sqrt(sf[i] * eta[i]**3 / (2 * tf[i] * rho * g))

print()
print(r)
print(eta)
print(q)
print(ms.mean_value(q))
print(ms.std_dev_m(q))
print(ms.std_dev_e(q))
