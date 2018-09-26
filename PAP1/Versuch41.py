import math as m
import measure as ms

# Gasthermometer
p0 = 1013.25
pl = 1020.2
dpl = 0.2

p = [916, 945, 991, 1027, 1061, 1092, 1125, 1163, 198, 1226, 1245]
dp = [1, 5, 6, 5, 5, 5, 3, 2, 5, 4, 2]
pc = 653
dpc = 2
pn = 255
dpn = 1


pnb = p[len(p) - 1] * p0 / pl
dpnb = p0 / pl * m.sqrt(dp[len(p) - 1]**2 + p[len(p) - 1]**2 / pl**2 * dpl**2)

ms.pve("pnb", pnb, dpnb)
