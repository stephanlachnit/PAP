import math as m
import measure as ms

# Messung mit variabler Länge
b0 = 42
db0 = 3
b = 30
db = 3
l = [1.5, 3.0, 6.0, 12.0, 24.0]

il = [[63280.29, 63070.40, 63058.95, 62993.61, 63105.31],
      [46607.18, 46647.68, 46587.64, 46623.69, 46628.55],
      [25813.76, 25864.76, 25832.21, 25804.32, 25801.63],
      [5206.06, 5182.66, 5183.84, 5176.01, 5202.05],
      [380.84, 389.33, 388.31, 387.34, 387.39]]
dil = [ms.std_dev_m(il[i]) for i in range(len(il))]
il = [ms.mean_value(il[i]) for i in range(len(il))]

il[4] *= (b * b) / (b0 * b0)
dil[4] = m.sqrt((4 * b / b0**2 * il[4] * db)**2 + (b**2 / b0**3 * il[4] * db0)**2 + (b**2 / b0**2 * dil[4])**2)

for i in range(len(l)):
  ms.pv("l" + str(i), l[i])
print()
for i in range(len(il)):
  ms.pve("il" + str(i), il[i], dil[i])

# Messung mit variaabler Konzentration
v = [21.0, 1.4, 1.6, 4.0, 14.0]
ct = 1.0e-3

ik = [[62905.49, 62315.77, 62290.07, 62279.29, 62243.04],
      [42832.73, 42892.47, 43043.32, 43101.21, 43228.01],
      [27760.81, 27836.17, 27755.93, 27715.35, 27603.68],
      [8208.01, 8366.15, 8491.14, 8561.95, 8646.58],
      [1487.15, 1458.05, 1564.75, 1531.03, 1604.25]]
dik = [ms.std_dev_m(ik[i]) for i in range(len(ik))]
ik = [ms.mean_value(ik[i]) for i in range(len(ik))]

c = []
vsum = 0
for i in range(len(ik)):
  c.append(ct * vsum / (vsum + v[0]))
  if i < len(ik) - 1:
    vsum += v[i + 1]

print()
print()
for i in range(len(c)):
  ms.pv("c" + str(i), c[i])
print()
for i in range(len(ik)):
  ms.pve("ik" + str(i), ik[i], dik[i])
