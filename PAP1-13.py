import measure as ms

# Constants
I_1 = 280e-3
I_2 = 360e-3
fg = 2500

# Measured values in SI-units
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

f2 = [302.99, 500.90, 702.15, 900.45, 1100.40, 1308.60, 1503.80, 1700.00, 1900.20, 2100.20, 1150.60, 1201.60, 1251.20, 1273.50, 1352.30, 1403.20]
df2 = [0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.10, 0.20, 0.20, 0.30, 0.20, 0.20, 0.20, 0.30, 0.20, 0.20]
phiE2 = [0.60, 0.60, 0.70, 1.00, 1.90, 3.60, 1.10, 0.60, 0.40, 0.20, 2.60, 3.70, 4.80, 4.60, 2.60, 1.80]
dPhiE2 = [0.10, 0.05, 0.10, 0.05, 0.10, 0.10, 0.10, 0.05, 0.10, 0.10, 0.10, 0.10, 0.05, 0.05, 0.10, 0.10]

# Graphical determinated values
n1 = 3.8
dn1 = 0.3
n2 = 2.4
dn2 = 0.2

fl1 = 1213
dfl1 = 13
fr1 = 1293
dfr1 = 13
fl2 = 1187
dfl2 = 13
fr2 = 1307
dfr2 = 13

bmax1 = 8.3
dbmax1 = 0.3
bmax2 = 4.9
dbmax2 = 0.1

# Measure series reduction
dT = ms.std_dev_m(T) / 20
T = ms.mean_value(T) / 20
dT1 = ms.std_dev_m(T1) / 15
T1 = ms.mean_value(T1) / 15
dT2 = ms.std_dev_m(T2) / 10
T2 = ms.mean_value(T2) / 10

print()
print(ms.val("T0", T, dT))
print(ms.val("T1", T1, dT1))
print(ms.val("T2", T2, dT2))

# Damping constant calculated with time of half amplitude
t1 = n1 * T1
dt1 = ms.sqrt((n1 * dT1)**2 + (T1 * dn1)**2)
t2 = n2 * T2
dt2 = ms.sqrt((n2 * dT2)**2 + (T2 * dn2)**2)
delta_ha1 = ms.ln(2) / t1
dDelta_ha1 = ms.ln(2) / t1**2 * dt1
delta_ha2 = ms.ln(2) / t2
dDelta_ha2 = ms.ln(2) / t2**2 * dt2

print()
print(ms.val("t1", t1, dt1))
print(ms.val("t2", t2, dt2))
print()
print(ms.val("detla_ha1", delta_ha1, dDelta_ha1))
print(ms.val("detla_ha2", delta_ha2, dDelta_ha2))

# Damping constant calculated with full width at half maximum
wl1 = 2 * ms.pi * fl1 / fg
dwl1 = 2 * ms.pi * dfl1 / fg
wr1 = 2 * ms.pi * fr1 / fg
dwr1 = 2 * ms.pi * dfr1 / fg
wl2 = 2 * ms.pi * fl2 / fg
dwl2 = 2 * ms.pi * dfl2 / fg
wr2 = 2 * ms.pi * fr2 / fg
dwr2 = 2 * ms.pi * dfr2 / fg

delta_hm1 = (wr1 - wl1) / 2
dDelta_hm1 = 1 / 2 * ms.sqrt(dwr1**2 + dwl1**2)
delta_hm2 = (wr2 - wl1) / 2
dDelta_hm2 = 1 / 2 * ms.sqrt(dwr2**2 + dwl2**2)

print()
print(ms.val("wl1", wl1, dwl1))
print(ms.val("wr1", wr1, dwr1))
print(ms.val("wl2", wl2, dwl2))
print(ms.val("wr2", wr2, dwr2))
print()
print(ms.val("delta_hm1", delta_hm1, dDelta_hm1))
print(ms.val("delta_hm2", delta_hm2, dDelta_hm2))

# Damping constant calculated with resonance magnification
[tmp, tmp, bn1, dbn1] = ms.linreg([f1[i] for i in range(3)], [phiE1[i] for i in range(3)], [dPhiE1[i] for i in range(3)], [df1[i] for i in range(3)])
[tmp, tmp, bn2, dbn2] = ms.linreg([f2[i] for i in range(3)], [phiE2[i] for i in range(3)], [dPhiE2[i] for i in range(3)], [df2[i] for i in range(3)])

w0 = 2 * ms.pi / T
dw0 = 2 * ms.pi / T**2 * dT

delta_rm1 = w0 * bn1 / (2 * bmax1)
dDelta_rm1 = 1 / (2 * bmax1) * ms.sqrt((w0 * dbn1)**2 + (bn1 * dw0)**2 + (w0 * bn1 * dbmax1 / bmax1)**2)
delta_rm2 = w0 * bn2 / (2 * bmax2)
dDelta_rm2 = 1 / (2 * bmax2) * ms.sqrt((w0 * dbn2)**2 + (bn2 * dw0)**2 + (w0 * bn2 * dbmax2 / bmax2)**2)

print()
print(ms.val("bn1", bn1, dbn1))
print(ms.val("bmax1", bmax1, dbmax1))
print(ms.val("bn2", bn2, dbn2))
print(ms.val("bmax2", bmax2, dbmax2))
print()
print(ms.val("delta_rm1", delta_rm1, dDelta_rm1))
print(ms.val("delta_rm2", delta_rm2, dDelta_rm2))

# Comparision
print("\n1:")
print(ms.sig("ha vs hm", delta_ha1, delta_hm1, dDelta_ha1, dDelta_hm1))
print(ms.sig("ha vs rm", delta_ha1, delta_rm1, dDelta_ha1, dDelta_rm1))
print(ms.sig("hm vs rm", delta_hm1, delta_rm1, dDelta_hm1, dDelta_rm1))
print(ms.val("mean value", 1 / 3 * (delta_ha1 + delta_hm1 + delta_rm1), 1 / 3 * ms.sqrt(dDelta_ha1**2 + dDelta_hm1**2 + dDelta_rm1**2)))

print("\n2:")
print(ms.sig("ha vs hm", delta_ha2, delta_hm2, dDelta_ha2, dDelta_hm2))
print(ms.sig("ha vs rm", delta_ha2, delta_rm2, dDelta_ha2, dDelta_rm2))
print(ms.sig("hm vs rm", delta_hm2, delta_rm2, dDelta_hm2, dDelta_rm2))
print(ms.val("mean value", 1 / 3 * (delta_ha2 + delta_hm2 + delta_rm2), 1 / 3 * ms.sqrt(dDelta_ha2**2 + dDelta_hm2**2 + dDelta_rm2**2)))
