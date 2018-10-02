import measure as ms

# Constants


# Measured values in SI-units
I_1 = 280e-3
I_2 = 360e-3
TE1 = 29.72
TE2 = 19.76

T = [39.78, 39.75, 39.78]

T1 = [29.72, 29.53, 29.72]
phiF1 = [20.0, 17.0, 14.2, 11.8, 9.8, 8.2, 7.0, 5.8, 4.8, 3.9, 3.2, 2.6, 2.2, 1.8, 1.4, 1.1, 0.9, 0.7]
dPhiF1 = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
T2 = [19.76, 19.75, 19.80]
phiF2 = [20.0, 14.8, 11.0, 8.1, 6.1, 4.6, 3.4, 2.6, 2.0, 1.5, 1.2, 1.0, 0.8, 0.6]
dPhiF2 = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]

f1 = [300.60, 501.18, 699.94, 901.19, 1101.00, 1200.40, 1301.40, 1500.80, 1700.50, 1907.20, 2101.30, 1151.00, 1250.90, 1351.00, 1401.30, 1271.10]
df1 = [0.10, 0.04, 0.05, 0.04, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.20, 0.20, 0.20]
phiE1 = [0.55, 0.60, 0.70, 1.00, 2.00, 4.60, 5.25, 1.10, 0.60, 0.40, 0.30, 2.80, 8.00, 2.30, 2.00, 7.60]
dPhiE1 = [0.05, 0.05, 0.10, 0.05, 0.05, 0.05, 0.10, 0.05, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10, 0.05, 0.10]
f2 =  [302.99, 500.90, 702.15, 900.45, 1100.40, 1308.60, 1503.80, 1700.00, 1900.20, 2100.20, 1150.60, 1201.60, 1251.20, 1273.50, 1352.30, 1403.20]
df2 = [0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.10, 0.20, 0.20, 0.30, 0.20, 0.20, 0.20, 0.30, 0.20, 0.20]
phiE2 = [0.60, 0.60, 0.70, 1.00, 1.90, 3.60, 1.10, 0.60, 0.40, 0.20, 2.60, 3.70, 4.80, 4.60, 2.60, 1.80]
dPhiE2 = [0.10, 0.05, 0.10, 0.05, 0.10, 0.10, 0.10, 0.05, 0.10, 0.10, 0.10, 0.10, 0.05, 0.05, 0.10, 0.10]

# Measure series reduction
dT = ms.std_dev_m(T) / 20
T = ms.mean_value(T) / 20

dT1 = ms.std_dev_m(T1) / 15
T1 = ms.mean_value(T1) / 15
dT2 = ms.std_dev_m(T2) / 10
T2 = ms.mean_value(T2) / 10

ms.ps("T", T, T1, dT, dT1)

#d = 2.0 * ms.pi * ms.sqrt((1.0 / T)**2 - (2.0 / T1)**2)

print()
ms.pve("T0", T, dT)
ms.pve("T1", T1, dT1, False)
ms.pve("T2", T2, dT2)

###
w0 = 2 * ms.pi * 1247 / 2500
dw0 = 2 * ms.pi * 13 / 2500
wf = [2 * ms.pi / T1 / 2500, 2 * ms.pi / T2 / 2500 ]
dwf = [2 * ms.pi / T1**2 / 2500, 2 * ms.pi / T2**2 / 2500]
[n1, n2] = [3.8, 2.4]
[dn1, dn2] = [0.3, 0.2]
delta_n = [ms.ln(2) / (n1 * T1), ms.ln(2) / (n2 * T2)]
dDelta_n = [ms.ln(2) / (n1 * T1) * ms.sqrt((dn1 / n1)**2 + (dT1 / T1)**2), ms.ln(2) / (n2 * T2) * ms.sqrt((dn2 / n2)**2 + (dT2 / T2)**2)]
delta_f = dDelta_f = [0.0, 0.0]
for i in range(len(wf)):
  delta_f[i] = ms.sqrt(w0**2 - wf[i]**2)
  dDelta_f[i] = ms.sqrt(((w0 * dw0)**2 + (wf[i] * dwf[i])**2) /(w0**2 - wf[i]**2))

#ms.ple("delta_n", delta_n, dDelta_n)
#ms.ple("delta_f", delta_f, dDelta_f)

###

#ms.pv("d", d)

# Damping constant calculation by amplitude measurement
n1 = [i for i in range(len(phiF1))]
lgPhi1 = [ms.lg(phiF1[i]) for i in range(len(phiF1))]
dLgPhi1 = [1 / ms.ln(10.0) * abs(dPhiF1[i] / phiF1[i]) for i in range(len(phiF1))]
n2 = [i for i in range(len(phiF2))]
lgPhi2 = [ms.lg(phiF2[i]) for i in range(len(phiF2))]
dLgPhi2 = [1 / ms.ln(10.0) * abs(dPhiF2[i] / phiF2[i]) for i in range(len(phiF2))]

[s1, ds1, b1, db1] = ms.linreg(n1, lgPhi1, dLgPhi1)
[s2, ds2, b2, db2] = ms.linreg(n2, lgPhi2, dLgPhi2)

d1 = - s1 / T1
dd1 = 1.0 / T1 * ms.sqrt(ds1**2 + (s1 / T1 * dT1)**2)
d2 = - s2 / T2
dd2 = 1.0 / T2 * ms.sqrt(ds2**2 + (s2 / T2 * dT2)**2)

#ms.pve("d1", d1, dd1, False)
#ms.pve("d2", d2, dd2)

nh1 = 3.90
dnh1 = 0.25
nh2 = 2.45
dnh2 = 0.05
delta11 = ms.ln(2.0) / (nh1 * T1)
dDelta11 = ms.ln(2.0) / (nh1 * T1) * ms.sqrt((dnh1 / nh1)**2 + (dT1 / T1)**2)
delta12 = ms.ln(2.0) / (nh2 * T2)
dDelta12 = ms.ln(2.0) / (nh2 * T2) * ms.sqrt((dnh2 / nh2)**2 + (dT2 / T2)**2)

#ms.pve("nh1", nh1, dnh1, False)
#ms.pve("nh2", nh2, dnh2, False)
ms.pve("delta11", delta11, dDelta11, False)
ms.pve("delta12", delta12, dDelta12, False)

# Calculation of characteristic quantities of the resonance curve
f0 = 1.0 / T
df0 = abs(dT / T**2)
omega0 = 2.0 * ms.pi * f0
dOmega0 = 2.0 * ms.pi * df0

fgMax1 = 1260.0
dfgMax1 = 5.0
hg1 = 85.0
dhg1 = 5.0
reg1 = 8.075 / 0.55
dreg1 = 1.0 / 0.55 * ms.sqrt(0.05**2 + (8.075 / 0.55 * 0.025)**2)

omegaMax1 = 2.0 * ms.pi * fgMax1 / 2500
dOmegaMax1 = 2.0 * ms.pi * dfgMax1 / 2500
h1 = 2.0 * ms.pi * hg1 / 2500
dh1 = 2.0 * ms.pi * dhg1 / 2500
re1 = reg1
dre1 = dreg1

delta21 = h1 / 2.0
dDelta21 = dh1 / 2.0

delta31 = omega0 / (2.0 * re1)
dDelta31 = 1.0 / (2.0 * re1) * ms.sqrt(dOmega0**2 + (omega0 * dre1 / re1)**2)

#ms.ps("omega0,omegaMax", omega0, omegaMax, dOmega0, dOmegaMax)
ms.pve("delta21", delta21, dDelta21, False)
ms.pve("delta31", delta31, dDelta31, False)

fgMax2 = 1255.0
dfgMax2 = 5.0
hg2 = 130.0
dhg2 = 5.0
reg2 = 4.825 / 0.6
dreg2 = 1.0 / 0.6 * ms.sqrt(0.025**2 + (4.825 / 0.6 * 0.025)**2)

omegaMax2 = 2.0 * ms.pi * fgMax2 / 2500
dOmegaMax2 = 2.0 * ms.pi * dfgMax2 / 2500
h2 = 2.0 * ms.pi * hg2 / 2500
dh2 = 2.0 * ms.pi * dhg2 / 2500
re2 = reg2
dre2 = dreg2

delta22 = h2 / 2.0
dDelta22 = dh2 / 2.0

delta32 = omega0 / (2.0 * re2)
dDelta32 = 1.0 / (2.0 * re2) * ms.sqrt(dOmega0**2 + (omega0 * dre2 / re2)**2)

ms.pve("delta22", delta22, dDelta22, False)
ms.pve("delta32", delta32, dDelta32)

# Plot of the resonance curve
#omega1 = [2.0 * ms.pi * f1[i] for i in range(len(f1))]
#dOmega1 = [2.0 * ms.pi * df1[i] for i in range(len(f1))]
#omega2 = [2.0 * ms.pi * f2[i] for i in range(len(f2))]
#dOmega2 = [2.0 * ms.pi * df2[i] for i in range(len(f2))]
#
#title = "Deflection phi in dependency of the angular frequency omega for I1"
#ms.plt.figure(2)
#ms.plot(title, "Omega", "Phi", omega1, phiE1, dPhiE1, dOmega1)
#
#title = "Deflection phi in dependency of the angular frequency omega for I2"
#ms.plt.figure(3)
#ms.plot(title, "Omega", "Phi", omega2, phiE2, dPhiE2, dOmega2)