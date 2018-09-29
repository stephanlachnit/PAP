import measure as ms

import numpy as np
import matplotlib.pyplot as plt

m = [0, 40, 80, 120, 160, 200, 240]

# Tatsächliche Messwerte
#dm = [11.5, 12.0, 10.0, 15.0, 14.4, 17.5, 12.6]
#phi = [0.0, 48.0, 105.0, 147.0, 184.0, 245.0, 303.0]
#dPhi = [10.0, 15.0, 17.0, 13.0, 13.0, 11.0, 13.0]

# Testwerte
dm = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
phi = [0.0, 48.0, 97.0, 147.0, 194.0, 245.0, 293.0]
dPhi = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

b = ms.reg_itc(m, phi, dm, dPhi)
db = ms.reg_itc_err(m, phi, dm, dPhi)
g = ms.reg_grad(m, phi, dm, dPhi, True)
dg = ms.reg_grad_err(m, phi, dm, dPhi)

n = 100
xBeg = 0
xEnd = 250
x = []
y0 = []
y1 = []
for i in range(n):
  x.append(xBeg + (xEnd - xBeg) / (n - 1) * i)
  y0.append( b + g * x[i])
  y1.append((b - db) + (g + dg) * x[i])

x = np.array(x)
y0 = np.array(y0)
y1 = np.array(y1)

m = np.array(m)
dm = np.array(dm)
phi = np.array(phi)
dPhi = np.array(dPhi)

plt.title("Tangential force inducing mass m - deflection angle phi - dependence\nof a rotating pendulum")
plt.xlabel("m")
plt.ylabel("phi")
plt.plot(x, y0, label="line of best fit")
plt.plot(x, y1, label="line of best fit error")
plt.errorbar(m, phi, yerr = dPhi, xerr = dm, fmt='x', capsize=2.0)
plt.legend()
plt.show()