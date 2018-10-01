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
T = [T[i] / 20.0 for i in range(len(T))]
dT = ms.std_dev_m(T)
T = ms.mean_value(T)

T1 = [T1[i] / 15 for i in range(len(T1))]
T2 = [T2[i] / 10 for i in range(len(T2))]
dT1 = ms.std_dev_m(T1)
T1 = ms.mean_value(T1)
dT2 = ms.std_dev_m(T2)
T2 = ms.mean_value(T2)

print()
ms.pve("T0", T, dT)
ms.pve("T1", T1, dT1, False)
ms.pve("T2", T2, dT2)

# Damping constant calculation by amplitude measurement
n1 = [i for i in range(len(phiF1))]
n2 = [i for i in range(len(phiF2))]
lnPhiF1 = [ms.ln(phiF1[i]) for i in range(len(phiF1))]
dLnPhiF1 = [abs(dPhiF1[i] / phiF1[i]) for i in range(len(phiF1))]
lnPhiF2 = [ms.ln(phiF2[i]) for i in range(len(phiF2))]
dLnPhiF2 = [abs(dPhiF2[i] / phiF2[i]) for i in range(len(phiF2))]

ms.plt.figure(1)
ms.plot_linreg("title", "xlabel", "ylabel", n2, [1e-20 for i in range(len(n2))], lnPhiF2, dLnPhiF2)

[s1, ds1, b1, db1] = ms.linreg(n1, lnPhiF1, dLnPhiF1)
[s2, ds2, b2, db2] = ms.linreg(n1, lnPhiF1, dLnPhiF1)

th1 = (ms.ln(phiF1[0] / 2.0) - b1) / s1
dth1 = 1.0 / abs(s1) * ms.sqrt((2.0 * dPhiF1[0] / phiF1[0])**2 + db1**2 + ((ms.ln(phiF1[0] / 2.0) - b1) * ds1 / s1)**2)
th2 = (ms.ln(phiF2[0] / 2.0) - b2) / s2
dth2 = 1.0 / abs(s2) * ms.sqrt((2.0 * dPhiF2[0] / phiF2[0])**2 + db2**2 + ((ms.ln(phiF2[0] / 2.0) - b2) * ds2 / s2)**2)
delta1 = ms.ln(2.0) / th1
dDelta1 = ms.ln(2.0) * abs(dth1 / th1**2)
delta2 = ms.ln(2.0) / th2
dDelta2 = ms.ln(2.0) * abs(dth2 / th2**2)

ms.pl("n1", n1)
ms.ple("phiF1", phiF1, dPhiF1)
ms.pl("n2", n2)
ms.ple("phiF2", phiF2, dPhiF2)
#ms.pve("th1", th1, dth1, False)
#ms.pve("th2", th2, dth2, False)
#ms.pve("delta1", delta1, dDelta1, False)
#ms.pve("delta2", delta2, dDelta2, False)

# Calculation of characteristic quantities of the resonance curve
f0 = 1.0 / T
df0 = abs(dT / T**2)
omega0 = 2.0 * ms.pi * f0
dOmega0 = 2.0 * ms.pi * df0

#fgMax = 
#dfgMax = 
#fMax = fgMax / 2500
#dfMax = dfgMax / 2500

#h = 
#dh = 
#re = 
#dre = 

#delta2 = h / 2.0
#dDelta2 = dh / 2.0

#delta3 = omega0 / (2.0 * re)
#dDelta3 = 1.0 / (2.0 * re) * ms.sqrt(dOmega0**2 + (omega0 * dre / re)**2)

#ms.ps("f0,fMax", f0, fMax, df0, dfMax)

ms.ple("f1", f1, df1)
ms.ple("phiE1", phiE1, dPhiE1)
ms.ple("f2", f2, df2)
ms.ple("phiE2", phiE2, dPhiE2)

# Plot of the resonance curve
omega1 = [2.0 * ms.pi * f1[i] for i in range(len(f1))]
dOmega1 = [2.0 * ms.pi * df1[i] for i in range(len(f1))]
omega2 = [2.0 * ms.pi * f2[i] for i in range(len(f2))]
dOmega2 = [2.0 * ms.pi * df2[i] for i in range(len(f2))]

title = "Deflection phi in dependency of the angular frequency omega for I1"
ms.plt.figure(2)
ms.plot(title, "Omega", "Phi", omega1, phiE1, dPhiE1, dOmega1)

title = "Deflection phi in dependency of the angular frequency omega for I2"
ms.plt.figure(3)
ms.plot(title, "Omega", "Phi", omega2, phiE2, dPhiE2, dOmega2)

# Das Problem ist, dass man die Werte sortieren muss. Viel Spaß Jona:-P
