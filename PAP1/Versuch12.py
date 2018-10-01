import measure as ms

zero = 1.0e-20

# Constants
g = 9.80984
dg = 2e-5

# Measured values in SI-units
r = 100e-3 / 2
dr = 0.5e-3 / 2
m = [0, 40e-3, 80e-3, 120e-3, 160e-3, 200e-3, 240e-3]
dm = [zero, 0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.6e-3]
phi = [0.0, 48.0, 97.0, 147.0, 194.0, 245.0, 293.0]
dPhi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

md = 552.3e-3
dmd = 0.1e-3
rd = 104.9e-3 / 2.0
drd = 0.5e-3 / 2.0
T1 = [24.10, 23.98, 24.18, 24.02, 24.07]
T2 = [33.15, 33.50, 33.23, 33.38, 33.09]

mbp = 679.00e-3
dmbp = 0.20e-3
T = [[46.56, 46.64, 46.62], [46.75, 46.73, 46.64], [47.54, 47.39, 47.39], [48.39, 48.20, 48.32], [49.96, 50.17, 50.09], [52.50, 52.48, 52.46]]
a = [0.0e-3, 5.7e-3, 11.8e-3, 18.6e-3, 29.9e-3, 36.8e-3]
da = [2.0e-3, 2.0e-3, 2.0e-3, 2.0e-3, 2.0e-3, 2.0e-3]

# Measurement series reduction
dT1 = ms.std_dev_m(T1) / 20
T1 = ms.mean_value(T1) / 20
dT2 = ms.std_dev_m(T2) / 20
T2 = ms.mean_value(T2) / 20
dT = [ms.std_dev_m(T[i]) / 20 for i in range(len(T))]
T = [ms.mean_value(T[i]) / 20 for i in range(len(T))]


# 1.1 Determination of the deflecting force by torque-deflection angle-dependency measurement
M = []
dM = []
for i in range(len(m)):
  M.append(m[i] * g * r)
  dM.append(ms.sqrt((m[i] * g * dr)**2 + (m[i] * r * dg)**2 + (g * r * dm[i])**2))

Dt = ms.reg_grad(phi, M, dPhi, dM) * 180.0 / ms.pi
dDt = ms.reg_grad_err(phi, M, dPhi, dM) * 180.0 / ms.pi

#ms.ple("ϕ", phi, dPhi)
#ms.ple("M", M, dM)
#ms.pve("D", Dt, dDt)
#print()

#ms.plot("Deflecting force determination", "Deflection angle ϕ", "Torque M", phi, dPhi, M, dM)


# 1.2 Deflecting force determination by period time measurement
Jd = md / 2.0 * rd**2
dJd = rd * ms.sqrt((0.5 * rd * dmd)**2 + (md * drd)**2)

Dp = 4.0 * ms.pi**2 * Jd / (T2**2 - T1**2)
dDp = 4.0 * ms.pi**2 / (T2**2 - T1**2) * ms.sqrt(dJd**2 + ((2.0 * Jd * T1 * dT1)**2 + (2.0 * Jd * T2 * dT2)**2) / (T2**2 - T1**2)**2)

#ms.pve("T1", T1, dT1, False)
#ms.pve("T2", T2, dT2, False)
#print(J)
#ms.pve("Jd", Jd, dJd, False)
#ms.pve("D", Dp, dDp)
#print()

# Deviation between 1.1 and 1.2
dD = max([dDt, dDp])
sigmaD = abs(Dt - Dp) / dD

ms.ps("D", Dt, Dp, dDt, dDp)


# 2 Verification of the theorem of Steiner by moment of inertia-center of mass distance-dependence measurement
a2 = []
da2 = []
J = []
dJ = []
Jt = []
dJt = []
for i in range(len(T)):
  a2.append(a[i]**2)
  da2.append(2.0 * abs(a[i] * da[i]))
  J.append(Dt / (4.0 * ms.pi**2) * (T[i]**2 - T1**2))
  dJ.append(ms.sqrt((Dt * T1 * dT1)**2 + (Dt * T[i] * dT[i])**2 + ((T[i]**2 - T1**2) * dDt)**2 / 4.0) / (2.0 * ms.pi**2))
  Jt.append(J[0] + mbp * a2[i])
  dJt.append(ms.sqrt(dJ[0]**2 + (a2[i] * dmbp)**2 + (mbp * da2[i])**2))

ms.ple("T", T, dT)
ms.ple("a2", a2, da2)
ms.ple("J", J, dJ)
ms.ple("dJt", Jt, dJt)

a2.pop(0)
da2.pop(0)
J.pop(0)
dJ.pop(0)
Jt.pop(0)
dJt.pop(0)

for i in range(len(a2)):
  ms.ps("J", J[i], Jt[i], dJ[i], dJt[i], False)
print()
for i in range(len(a2)):
  ms.pv("p", abs(1.0 - J[i] / Jt[i]) * 100.0, False)
print()

chi2 = ms.chi2(a2, da2, Jt, J, [max(dJ[i], dJt[i]) for i in range(len(a2))])
ms.pv("chi2", chi2)