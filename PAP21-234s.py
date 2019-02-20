# measure version: 1.3
from measure import h,e,c,val,lst,sqrt,plot,sig,chi2,chi2_red

# Frauenhoferlinien
# Measured Values
wl = [393e-9, 397e-9, 430e-9, 486e-9, 517e-9, 527e-9, 656e-9, 687e-9]
dwl = 1e-9
fwl = [393.4e-9, 396.8e-9, 430.8e-9, 486.1e-9, 518.4e-9, 527.0e-9, 656.3e-9, 686.7e-9]

print()
print("Frauenhoferlinien:")
for i in range(len(wl)):
  print(sig("{0:03.1f}".format(fwl[i]*1e9)+"nm", wl[i], dwl, fwl[i]))

print()
print("Nicht eindeutig:")
print(sig("587.6nm", 589e-9, dwl, 587.6e-9))
print(sig("589.0nm", 589e-9, dwl, 589.0e-9))
print(sig("589.6nm", 589e-9, dwl, 589.6e-9))

# Natrium
# Measured Values
wl = [331e-9, 395e-9, 404e-9, 416e-9, 420e-9, 427e-9, 430e-9, 434e-9, 439e-9, 442e-9, 449e-9, 454e-9, 467e-9, 475e-9, 498e-9, 515e-9, 525e-9, 568e-9, 589e-9, 668e-9, 675e-9, 687e-9, 696e-9, 703e-9, 707e-9, 715e-9, 720e-9, 727e-9, 738e-9, 751e-9, 763e-9, 770e-9, 773e-9, 795e-9, 801e-9, 811e-9, 819e-9]
dwl = [1.5e-9 for i in wl]

# Wellenlängen 1. Nebenserie
ns1 = [819e-9]
m = 3
Eryd = -13.605
E3p = Eryd / m**2 - h / e * c / ns1[0]

for m in range(4,13):
  ns1.append(h / e * c / (Eryd * m**-2 - E3p))

print()
print(lst("Wellenlängen 1. Nebenserie", ns1, [0.3e-9 for i in range(len(ns1))]))

# Wellenlängen 2. Nebenserie
ns2 = 589e-9
E3s = E3p - h / e * c / ns2
sCorr = 3 - sqrt(Eryd / E3s)

ns2 = []
for m in range(5,10):
  ns2.append(h / e * c / (Eryd  * (m - sCorr)**-2 - E3p))

print()
print(lst("Wellenlängen 2. Nebenserie", ns2, [0.3e-9 for i in range(len(ns2))]))

# Hauptserie
pCorr = 3 - sqrt(Eryd / E3p)

hs = []
for m in range(4,6):
  hs.append(h / e * c / (Eryd  * (m - pCorr)**-2 - E3s))

print()
print(lst("Wellenlängen Hauptserie", hs, [0.3e-9 for i in range(len(hs))]))

# Plot
wlplot = plot(title="Wellenlängen", xlabel="Wellenlänge", ylabel="Wellenlänge")
wlplot.plotdata(ns1, ns1, label="1. Nebenserie")
wlplot.plotdata(ns2, ns2, label="2. Nebenserie")
wlplot.plotdata(hs, hs, label="Hauptserie")
wlplot.plotdata(wl, wl, dx=dwl, dy=dwl, label="Gemessen")
wlplot.drawplot()

# Identifizierungen
_dwl = dwl[0]
print()
print("Abweichungen:")
print(sig("HS  332.1nm", wl[0], _dwl, hs[0]))
print(sig("NS1 423.0nm", wl[4], _dwl, ns1[9]))
print(sig("NS1 425.6nm", wl[5], _dwl, ns1[8]))
print(sig("NS1 429.1nm", wl[6], _dwl, ns1[7]))
print(sig("NS1 433.9nm", wl[7], _dwl, ns1[6]))
print(sig("NS1 440.8nm", wl[8], _dwl, ns1[5]))
print(sig("NS2 444.1nm", wl[9], _dwl, ns2[4]))
print(sig("NS1 451.2nm", wl[10], _dwl, ns1[4]))
print(sig("NS2 456.5nm", wl[11], _dwl, ns2[3]))
print(sig("NS1 468.3nm", wl[12], _dwl, ns1[3]))
print(sig("NS2 477.6nm", wl[13], _dwl, ns2[2]))
print(sig("NS1 499.7nm", wl[14], _dwl, ns1[2]))
print(sig("NS2 518.7nm", wl[15], _dwl, ns2[1]))
print(sig("NS1 570.0nm", wl[17], _dwl, ns1[1]))
print(sig("NS1 819.0nm", wl[36], _dwl, ns1[0]))

# Bestimmung der Serienenergien und Korrekturfaktoren
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as sci_chi2

# 1. Nebenserie
print()
print("Fit für 1. Nebenserie:")
wl1 = [wl[36], wl[17], wl[14], wl[12], wl[10], wl[8], wl[7], wl[6], wl[5], wl[4]]
dwl1 = [_dwl for i in range(len(wl1))]
qz = [float(i) for i in range(3,13)]
def fit_func1(m,E_Ry,E_3p,D_d):
  return h / e * c / (E_Ry * (m - D_d)**-2 - E_3p)
para = [Eryd, E3p, pCorr]
popt, pcov = curve_fit(fit_func1, qz, wl1, sigma=dwl1, p0=para)
Ery1 = popt[0]
dEry1 = sqrt(pcov[0][0])
E3p1 = popt[1]
dE3p1 = sqrt(pcov[1][1])
print(val("Ery", Ery1, dEry1))
print(val("E3p", E3p1, dE3p1))
print(val("dCorr", popt[2], sqrt(pcov[2][2])))
_chi2 = chi2(wl1, dwl1, fit_func1(qz,*popt))
dof = len(qz) - 3
_chi2_red = chi2_red(wl1, dwl1, fit_func1(qz,*popt), dof=dof)
print(val("chi2", _chi2))
print(val("chi2_red", _chi2_red))
print(val("Wkeit", 1 - sci_chi2.cdf(_chi2, dof)))

wl1plot = plot(title="1. Nebenserie", xlabel="Quantenzahl", ylabel="Wellenlänge", figure=2)
wl1plot.plotdata(qz, wl1, dwl1, label="Messwerte")
x = np.linspace(2.8, 12.2, 100)
wl1plot.plotfunc(x, fit_func1(x,*popt), label="Fitfunktion")
wl1plot.drawplot()
wl1plot.saveplot("1. Nebenserie.pdf")

# 2. Nebenserie
print()
print("Fit für 2. Nebenserie:")
wl2 = [wl[15], wl[13], wl[11], wl[9]]
dwl2 = [_dwl for i in range(len(wl2))]
qz = [float(i) for i in range(6,10)]
def fit_func2(m,E_Ry,E_3p,D_s):
  return h / e * c / (E_Ry * (m - D_s)**-2 - E_3p)
para = [Eryd, E3p, sCorr]
popt, pcov = curve_fit(fit_func2, qz, wl2, sigma=dwl2, p0=para)
Ery2 = popt[0]
dEry2 = sqrt(pcov[0][0])
E3p2 = popt[1]
dE3p2 = sqrt(pcov[1][1])
print(val("Ery", Ery2, dEry2))
print(val("E3p", E3p2, dE3p2))
print(val("sCorr", popt[2], sqrt(pcov[2][2])))
_chi2 = chi2(wl2, dwl2, fit_func2(qz,*popt))
dof = len(qz) - 3
_chi2_red = chi2_red(wl2, dwl2, fit_func2(qz,*popt), dof=dof)
print(val("chi2", _chi2))
print(val("chi2_red", _chi2_red))
print(val("Wkeit", 1 - sci_chi2.cdf(_chi2, dof)))

wl2plot = plot(title="2. Nebenserie", xlabel="Quantenzahl", ylabel="Wellenlänge", figure=3)
wl2plot.plotdata(qz, wl2, dwl2, label="Messwerte")
x = np.linspace(3.8, 10.2, 100)
wl2plot.plotfunc(x, fit_func2(x,*popt), label="Fitfunktion")
wl2plot.drawplot()
wl2plot.saveplot("2. Nebenserie.pdf")

# Abweichung
print()
print("Abweichung der Fitwerte:")
print(sig("Ery", Ery1, dEry1, Ery2, dEry2))
print(sig("E3p", E3p1, dE3p1, E3p2, dE3p2))

plot.showplots()
