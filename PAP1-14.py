import measure as ms

# Constants
gl = 9.80984
dgl = 2e-5
rhoe = 7.86e3
rholn = 1.20
rf = 0.2e-3 / 2
pln = 101325.0

# Measured values
t0 = [39.39, 39.36, 39.31]
n0 = 20
lpo = [0.9900, 0.9900, 0.9950, 0.9810, 0.9850, 0.9860]
lpu = [0.0150, 0.0150, 0.0210, 0.0055, 0.0110, 0.0105]
dlp = 0.0005
rk = 30.0e-3 / 2.0
drk = 0.1e-3 / 2.0

pl = 100980.0
dpl = 60.0

t = 982.89
n = 500

s = [8e-3, 16e-3, 21e-3, 26e-3, 34e-3, 38e-3, 42e-3, 48e-3, 50e-3, 54e-3, 58e-3, 60e-3, 64e-3, 67e-3, 70e-3, 73e-3, 76e-3, 78e-3, 81e-3, 84e-3, 85e-3]
ds = [3e-3 for i in range(len(s))]
sn = 119e-3
dsn = 1e-3
sa = 852e-3
dsa = 3e-3

# Preparation
lp = [lpo[i] - lpu[i] for i in range(len(lpo))]
dlp = ms.std_dev_m(lp) # "+" ms.sqrt(2) * dlp # If one uses the addition of the systematic and statistic uncertainty, it must be quadratic
lp = ms.mean_value(lp)
l = lp - rk
dl = ms.sqrt(dlp**2 + drk**2)
dt = ms.std_dev_e(t0)
T0 = ms.mean_value(t0) / n0
n1 = 2 * l * dt / (0.3 * T0 * dl)
print()
print(ms.val("n1", int(n1 + 0.5)))
print()
print(ms.val("gl", gl, dgl))

# Simple calculation
T = t / n
dT = dt / n
gs = 4 * ms.pi**2 * l / T**2
dgs = 4 * ms.pi**2 / T**2 * ms.sqrt(dl**2 + (2 * l * dT / T)**2)

print()
print(ms.val("gs", gs, dgs))
print(ms.sig("dev gs & gl", gs, dgs, gl, dgl))

# Complex calculation
mk = rhoe * 4 / 3 * ms.pi * rk**3
dmk = rhoe * 4 * ms.pi * rk**2 * drk # Wrong uncertainty calculation
mf = rhoe * ms.pi * rf**2 * (lp - 2 * rk)
dmf = rhoe * ms.pi * rf**2 * ms.sqrt(dlp**2 + 4 * drk**2)

ts = [50 * T * i for i in range(len(s))]
dts = [50 * i * dT for i in range(len(ts))] # Wrong uncertainty calculation
s = [sn - s[i] for i in range(len(s))]
ds = [ms.sqrt(dsn**2 + ds[i]**2) for i in range(len(ds))]
phi = [ms.arctan(s[i] / sa) for i in range(len(s))]
dphi = [1 / (s[i]**2 + sa**2) * ms.sqrt((s[i] * dsa)**2 + (sa * ds[i])**2) for i in range(len(s))]
phi0 = ms.mean_value(phi)
dphi0 = ms.std_dev_m(phi) # I do not think one should calculate the uncertainty like that, because the deviation of the values is naturally very high
lnphi = [ms.ln(phi[i]) for i in range(len(s))]
dlnphi = [abs(dphi[i] / phi[i]) for i in range(len(s))] # The uncertainty here is just to be calculated with gaussian uncertainty propagation
[delta, dDelta, tmp, tmp] = ms.linreg(ts, lnphi, dlnphi, dts)
delta *= -1
ms.plot_linreg(ts, lnphi, dlnphi, dts)

omega0 = 2 * ms.pi / T
dOmega0 = 2 * ms.pi / T**2 * dT

rhol = pl / pln * rholn
drhol = dpl / pln * rholn
drhol = ms.sqrt(drhol**2 + ((25.0 + 273.15) / (20.0 + 273.15) * rhol - rhol)**2)

gc = 4 * ms.pi**2 * l / T**2 * (1 + 2 / 5 * (rk / l)**2 + rhol / rhoe - mf / (6 * mk) + (delta / omega0)**2 + phi0**2 / 8)
dgc = 4.0 * ms.pi**2 * l / T**2 * ms.sqrt((1.0 + 2.0 / 5.0 * (rk / l)**2 + rhol / rhoe - mf / (6.0 * mk) + (delta / omega0)**2 + phi0**2 / 8.0)**2 * ((dl / l)**2 + (2.0 * dT / T)**2)
  + (4.0 / 5.0 * (rk / l)**2)**2 * ((drk / rk)**2 + (dl / l)**2)
  + (drhol / rhoe)**2
  + (1.0 / 6.0 * mf / mk)**2 * ((dmf / mf)**2 + (dmk / mk)**2)
  + (2.0 * (delta / omega0)**2)**2 * ((dDelta / delta)**2 + (dOmega0 / omega0)**2)
  + (phi0 * dphi0 / 4)**2) # Wrong uncertainty calculation. Corrected calculation is formatted differently.

print()
print(ms.val("gc", gc, dgc))
print(ms.sig("dev gc & gl", gc, dgc, gl, dgl))

# Comparison
print()
print(ms.sig("dev gs & gc", gs, dgs, gc, dgc))
