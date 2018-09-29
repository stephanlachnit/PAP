import measure as ms

# constants
g = 9.80984
dg = 2e-5

# measured values in SI-units
m = [0e-3, 40e-3, 80e-3, 120e-3, 160e-3, 200e-3, 240e-3]
dm = [0.0e-3, 0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.6e-3]
phi = [0.0, 48.0, 97.0, 147.0, 194.0, 245.0, 293.0]
dPhi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
r = 100e-3 / 2.0
dr = 0.5e-3 / 2.0
T1 = [24.10, 23.98, 24.18, 24.02, 24.07]
T2 = [33.15, 33.50, 33.23, 33.38, 33.09]
rd = 104.9e-3 / 2.0
drd = 0.5e-3 / 2.0
md = 552.30e-3
dmd = 0.10e-3
T = [[46.56, 46.64, 46.62], [46.75, 46.73, 46.64], [47.54, 47.39, 47.39], [48.39, 48.20, 48.32], [49.96, 50.17, 50.09], [52.50, 52.48, 52.46]]
a = [0.0e-3, 5.7e-3, 11.8e-3, 18.6e-3, 29.9e-3, 36.8e-3]
da = [0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3, 0.6e-3]
mbp = 679.00e-3
dmbp = 0.20e-3

# measurement series reduction
dT1 = ms.std_dev_m(T1)
T1 = ms.mean_value(T1)
dT2 = ms.std_dev_m(T2)
T2 = ms.mean_value(T2)
dT = [ms.std_dev_m(T[i]) for i in range(len(T))]
T = [ms.mean_value(T[i]) for i in range(len(T))]

# calculation
M = []
dM = []
for i in range(len(m)):
  M.append(m[i] * g * r)
  dM.append(ms.sqrt((m[i] * g * dr)**2 + (m[i] * r * dg)**2 + (g * r * dm[i])**2))

tmp = ms.plot("Versuch 12", "Winkel ϕ", "Drehmoment M", phi, dPhi, M, dM)
D = tmp[0] * 180 / ms.pi()
dD = tmp [1] * 180 / ms.pi()

ms.ple("Winkel", phi, dPhi)
ms.ple("Drehmoment", M, dM)
ms.pve("Richtmoment", D, dD)

# 