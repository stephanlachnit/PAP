import measure as ms

# Constants
gl = 9.80984
dgl = 2e-5
rohe = 7.86e3
rohl = 1.20
rf = 0.2e-3 / 2

# Measured values
t0 = [39.39, 39.36, 39.31]
n0 = 20
lpo = [0.9900, 0.9900, 0.9950, 0.9810, 0.9850, 0.9860]
lpu = [0.0150, 0.0150, 0.0210, 0.0055, 0.0110, 0.0105]
dlp = 0.0005
rk = 30.0e-3 / 2.0
drk = 0.1e-3 / 2.0

pl = 1009.8e2
dpl = 0.6e2

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
dlp = ms.std_dev_m(lp) + ms.sqrt(2) * dlp
lp = ms.mean_value(lp)
l = lp - rk
dl = ms.sqrt(dlp**2 + drk**2)
dt = ms.std_dev_e(t0)
T0 = ms.mean_value(t0) / n0
n1 = 2 * l * dt / (0.3 * T0 * dl)
print()
print(ms.val("n1", int(n1 + 0.5)))

# Simple calculation
T = t / n
dT = dt / n
gs = 4 * ms.pi**2 * l / T**2
dgs = 4 * ms.pi**2 / T**2 * ms.sqrt(dl**2 + (2 * l * dT / T)**2)

print()
print(ms.val("T", T, dT))
print()
print(ms.val("g simple", gs, dgs))
print(ms.sig("dev gs & gl", gs, dgs, gl, dgl))

# Complex calculation
mk = rohe * 4 / 3 * ms.pi * rk**3
dmk = rohe * 4 * ms.pi * drk**2
mf = rohe * ms.pi * rf**2 * (lp - 2 * rk)
dmf = rohe * ms.pi * rf**2 * ms.sqrt(dlp**2 + 4 * drk**2)

ts = [50 * T * i for i in range(len(s))]
dts = [ts[i] * dT for i in range(len(ts))]
s = [sn - s[i] for i in range(len(s))]
ds = [ms.sqrt(dsn**2 + ds[i]**2) for i in range(len(ds))]
phi = [ms.arctan(s[i] / sa) for i in range(len(s))]
dphi = [1 / (s[i]**2 + sa**2) * ms.sqrt((s[i] * dsa)**2 + (sa * ds[i])**2) for i in range(len(s))]
phi0 = ms.mean_value(phi)
dphi0 = ms.std_dev_m(phi)
lnphi = [ms.ln(phi[i]) for i in range(len(s))]
lndphi = [ms.ln(phi[i] + dphi[i]) - ms.ln(phi[i]) for i in range(len(s))]
[delta, dDelta, tmp, tmp] = ms.linreg(ts, lnphi, lndphi, dts)
delta *= - 1

w0 = 2 * ms.pi / T
dw0 = 2 * ms.pi / T**2 * dT

gc = 4 * ms.pi**2 * l / T**2 * (1 + 2 / 5 * (rk / l)**2 + rohl / rohe - mf / (6 * mk) + (delta / w0)**2 + phi0**2 / 8)
dgc = 4 * ms.pi**2 / T**2 * ms.sqrt(l**2 * (4* delta**2 / w0**4 * ((delta * dw0 / w0)**2 + (l * dDelta)**2)
                                           + 1 / (36 * mk**2) * (dmf**2 + (mf * dmk / mk)**2)
                                           + 1 / 16 * (phi0 * dphi0)**2
                                           + ((1 + 2 / 5 * (rk / l)**2 + rohl / rohe - mf / (6 * mk) + (delta / w0)**2 + phi0**2 / 8) * dl)**2
                                           + ((1 - 2 / 5 * (rk / l)**2 + rohl / rohe - mf / (6 * mk) + (delta / w0)**2 + phi0**2 / 8) * 2 * dT)**2 )
                                    + (4 / 5 * rk * drk / l)**2)

print()
print(ms.val("mk", mk, dmk))
print(ms.val("mf", mf, dmf))
print()
print(ms.val("delta", delta, dDelta))
print(ms.val("phi0", phi0 * 180 / ms.pi, dphi0 * 180 / ms.pi))
print()
print(ms.val("g complex", gc, dgc))
print(ms.sig("dev gc & gl", gc, dgc, gl, dgl))

# Comparison
print()
print(ms.sig("dev gs & gc", gs, dgs, gc, dgc))
