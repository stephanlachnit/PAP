import measure as ms

# Richtmoment
g = 9.80984
dg = 2e-5
r = 100e-3 / 2
dr = 0.5e-3 / 2
m = [0, 40, 80, 120, 160, 200, 240]
dm = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
phi = [0.0, 48.0, 97.0, 147.0, 194.0, 245.0, 293.0]
dPhi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
M = []
dM = []
for i in range(len(m)):
  M.append(m[i] * g * r)
  dM.append(ms.sqrt((m[i] * g * dr)**2 + (m[i] * r * dg)**2 + (g * r * dm[i])**2))

ms.ple("Drehmoment", M, dM)
ms.plot("Versuch 42", "Drehmoment M", "Winkel φ", M, dM, phi, dPhi)
