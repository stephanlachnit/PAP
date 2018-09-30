import measure as ms

# Richtmoment
g = 9.80984
dg = 2e-5
r = 100e-3 / 2
dr = 0.5e-3 / 2
m = [0, 40e-3, 80e-3, 120e-3, 160e-3, 200e-3, 240e-3]
dm = [0.01e-3, 0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.6e-3]
phi = [0.0, 48.0, 97.0, 147.0, 194.0, 245.0, 293.0]
dPhi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
M = []
dM = []
for i in range(len(m)):
  M.append(m[i] * g * r)
  dM.append(ms.sqrt((m[i] * g * dr)**2 + (m[i] * r * dg)**2 + (g * r * dm[i])**2))

ms.plot("Deflecting force determination", "Deflection angle ϕ", "Torque M", phi, dPhi, M, dM)
D = ms.reg_grad(m, phi, dm, dPhi) * 180 / ms.pi()
dD = ms.reg_grad_err(m, phi, dm, dPhi) * 180 / ms.pi()

ms.ple("ϕ", phi, dPhi)
ms.ple("M", M, dM)
ms.pve("D", D, dD)

# Trägheitsmoment
m = 552.3e-3
dm = 0.1e-3
r = 104.9e-3 / 2.0
dr = 0.5e-3 / 2.0
J = m / 2.0 * r**2
dJ = r**2 / 2.0 * ms.sqrt((r * dm)**2 + (2.0 * m * dr)**2)

T1 = [24.10, 23.98, 24.18, 24.02, 24.07]
T2 = [33.15, 33.50, 33.23, 33.38, 33.09]
dT1 = ms.std_dev_m(T1) / 20
T1 = ms.mean_value(T1) / 20
dT2 = ms.std_dev_m(T2) / 20
T2 = ms.mean_value(T2) / 20
D = 2 * ms.pi()**2 * m * r**2 / (T2**2 - T1**2)

ms.pve("Trägheitsmoment", J, dJ)
ms.pve("Richtmoment", D, dD)
