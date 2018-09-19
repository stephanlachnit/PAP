import math as m
import measure as ms

# Abbildungsmaßstab
G = [15.0 for i in range(7)]
d_G = [1.0 for i in range(7)]
g = [210.0, 230.0, 240.0, 186.0, 170.0, 160.0, 150.0]
d_g = [1.0 for i in range(7)]
B = [13.5, 11.0, 10.0, 15.0, 22.0, 26.0, 28.0]
d_B = [1.0 for i in range(7)]
b = [190.0, 173.0, 166.0, 198.0, 240.0, 265.0, 280.0]
d_b = [10.0, 7.0, 7.0, 7.0, 20.0, 20.0, 20.0]
beta = [0.0 for i in range(7)]
d_beta = [0.0 for i in range(7)]
for i in range(7):
 beta[i] = ms.mean_value([B[i]/G[i], b[i]/g[i]])
 d_beta[i] = 0.5 * m.sqrt((d_B[i]**2 + (B[i] * d_G[i] / G[i])**2) / G[i]**2 + (d_b[i]**2 + (b[i] * d_g[i] / g[i])**2) / g[i]**2 )

#print()
#print(beta)
#print(d_beta)

# Besselverfahren
L = [0.7, 0.65, 0.6]
d_L = 1.0e-3
s1 = [0.352, 0.359, 0.366]
d_s1 = 1.0e-2
s2 = [0.747, 0.691, 0.635]
d_s2 = 5.0e-3
d = [0.0, 0.0, 0.0]
d_d = m.sqrt(d_s2**2 + d_s1**2)
f = [0.0, 0.0, 0.0]
d_f = [0.0, 0.0, 0.0]
for i in range(3):
 d[i] = s2[i] - s1[i]
 f[i] = (L[i]**2 - d[i]**2) / (4 * L[i])
 d_f[i] = m.sqrt(((L[i]**2 + d[i]**2) * d_L / (2 * L[i]))**2 + (d[i] * d_d)**2) / (2 * L[i])
f_mv = ms.mean_value(f)
d_f_mv = 1 / 3 * m.sqrt(d_f[0]**2 + d_f[1]**2 + d_f[2]**2)

#print()
#print(d)
#print(d_d)
#print()
#print(f)
#print(d_f)
#print()
#print(f_mv)
#print(d_f_mv)

# chromatische Aberration
L = 0.65
d_L = 1.0e-3
s1 = [[0.360, 0.363, 0.358], [0.353, 0.359, 0.35]]
d_s1 = [0.0, 0.0]
s2 = [[0.691, 0.69, 0.69], [0.693, 0.688, 0.693]]
d_s2 = [0.0, 0.0]
d = [0.0, 0.0]
d_d = [0.0, 0.0]
f = [0.0, 0.0]
d_f = [0.0, 0.0]
for i in range(2):
 d_s1[i] = ms.std_dev_m(s1[i])
 s1[i] = ms.mean_value(s1[i])
 d_s2[i] = ms.std_dev_m(s2[i])
 s2[i] = ms.mean_value(s2[i])
 d[i] = s2[i] - s1[i]
 d_d[i] = m.sqrt(d_s2[i]**2 + d_s1[i]**2)
 f[i] = (L**2 - d[i]**2) / (4 * L)
 d_f[i] = m.sqrt(((L**2 + d[i]**2) * d_L / (2 * L))**2 + (d[i] * d_d[i])**2) / (2 * L)
sigma = [abs(f[0] - f[1]) / d_f[0], abs(f[0] - f[1]) / d_f[1]]

#print()
#print(d)
#print(d_d)
#print()
#print(f)
#print(d_f)
#print(sigma)

# Auflösungsvermögen Mikroskop
s = 41.5
ds = 1.0
b = 250
f = 40
lam = 550e-6
beta = b / f - 1.0

a = [0.2, 0.5, 0.4]
da = ms.std_dev_m(a)
a = ms.mean_value(a)

alpha = m.atan(a / (2 * s))
dalpha = 1.0 / m.sqrt(4 * s**2 + a**2) * m.sqrt(da**2 + (a / s)**2 * ds**2)

g = 1 / (2 * beta)

gMin = 0.61 * lam / m.sin(alpha)
dgMin = abs(0.61 * lam * (m.cos(alpha) / (m.sin(alpha)**2)) * dalpha)

ms.pve("a", a, da)
ms.pv("beta", beta)
ms.pve("alpha", alpha * 180 / m.pi, dalpha * 180 / m.pi)
ms.pv("g", g)
ms.pve("gMin", gMin, dgMin)
ms.pv("sigma", abs(gMin - g) / dgMin)
