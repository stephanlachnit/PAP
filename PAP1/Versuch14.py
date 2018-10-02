import measure as ms

# Constants
gl = 9.80984
dgl = 2e-5
df = 7.86e3

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
T = t / n

s = [8e-3, 16e-3, 21e-3, 26e-3, 34e-3, 38e-3, 42e-3, 48e-3, 50e-3, 54e-3, 58e-3, 60e-3, 64e-3, 67e-3, 70e-3, 73e-3, 76e-3, 78e-3, 81e-3, 84e-3, 85e-3]
ds = [3e-3 for i in range(len(s))]
sn = 119e-3
dsn = 1e-3
sa = 852e-3
dsa = 3e-3
ts = [50 * T * i for i in range(len(s))]

# preparation
lp = [lpo[i] - lpu[i] for i in range(len(lpo))]
ms.pv("dlp", ms.sqrt(2) * dlp)
dlp = ms.std_dev_m(lp)
ms.pv("dlp", dlp)
lp = ms.mean_value(lp)
l = lp - rk
dl = ms.sqrt(dlp**2 + drk**2)
dt = ms.std_dev_e(t0)
T0 = ms.mean_value(t0) / n0
n1 = 2 * dt * l / (0.3 * T0 * dl)
ms.pv("n1", n1)

# period time
g = 4 * ms.pi**2 * l / T**2
dg = 4 * ms.pi**2 / T**2 * ms.sqrt(dl**2 + (2 * l * dt / t)**2)

ms.pve("g", g, dg)
ms.ps("dev g", gl, dgl, g, dg)
