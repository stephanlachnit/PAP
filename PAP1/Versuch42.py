import math as m
import measure as ms

p0 = 1013.0
cw = 4.186e3
dcw = 0.004e3

# Wasserwert
mk = 251.0e-3
dmk = 0.02e-3
mw = 497.61e-3
dmw = 0.13e-3
mw -= mk
dmw = m.sqrt(dmw**2 + dmk**2)
Tw = 55.8
dTw = 0.2
Tl = 24.1
dTl = 0.1
Tm = 53.25
dTm = 0.27
W = mw * cw * (Tw - Tm) / (Tm - Tl)
dW = m.sqrt((Tw - Tm)**2 * ((cw * dmw)**2 + (mw * dcw)**2) + (mw * cw)**2 * (dTw**2 + ( ((Tl - Tw) * dTm)**2 + ((Tw - Tm) * dTl)**2) / (Tm - Tl)**2)) / (Tm - Tl)

W = 92.0
dW = 17.0

ms.pve("Wasserwert", W, dW)

# T=100°C
p = 1024.6
dp = 2.0
mp = [501.773e-3, 131.08e-3, 125.05e-3]
dmp = [0.06e-3, 0.04e-3, 0.02e-3]
mw = [612.38e-3, 600.03e-3, 615.13e-3]
dmw = [0.04e-3, 0.09e-3, 0.05e-3]
Tw = [27.4, 24.7, 23.4]
dTw = [0.2, 0.2, 0.2]
Tp = [100 + 0.0276 * (p - p0) for i in range(len(mp))]
dTp = [0.0276 * dp for i in range(len(mp))]
Tm = [29.7, 29.5, 27.4]
dTm = [0.1, 0.1, 0.1]
cp = []
dcp = []
for i in range(len(mp)):
  mw[i] -= mk
  dmw[i] = m.sqrt(dmw[i]**2 + dmk**2)
  cp.append((mw[i] * cw + W) * (Tm[i] - Tw[i]) / (mp[i] * (Tp[i] - Tm[i])))
  dcp.append(1.0 / (mp[i] * (Tp[i] - Tm[i])) * m.sqrt(
      (Tm[i] - Tw[i])**2 * ((cw * dmw[i])**2 + (mw[i] * dcw)**2 + dW**2
    + (mw[i] * cw + W)**2 * ((dmp[i] / mp[i])**2 + (dTp[i] / (Tp[i] - Tm[i]))**2))
    + (mw[i] * cw + W)**2 * (dTw[i]**2 + ((Tp[i] - Tw[i]) / (Tp[i]- Tm[i]) * dTm[i])**2)))

sigma = [129.0, 900.0, 709.0]
for i in range(len(sigma)):
  sigma[i] = (sigma[i] - cp[i]) / dcp[i]

ms.ple("cp", cp, dcp)
ms.pl("sigma", sigma)

# T=-196°C
