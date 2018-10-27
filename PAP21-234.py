from measure import h,e,c,lst,sqrt,initplot,plot,showplot,sig

# Frauenhoferlinien
# Measured Values
wl = [393e-9, 397e-9, 430e-9, 486e-9, 517e-9, 527e-9, 656e-9, 687e-9]
dwl = 1.5e-9
fwl = [393.4e-9, 396.8e-9, 430.8e-9, 486.1e-9, 518.4e-9, 527.0e-9, 656.3e-9, 686.7e-9]

print()
print("Frauenhoferlinien")
print()
for i in range(len(wl)):
  print(sig("{0:03.1f}".format(fwl[i]*1e9)+"nm", wl[i], dwl, fwl[i]))

print(sig("587.6nm", 589e-9, dwl, 587.6e-9))
print(sig("589.0nm", 589e-9, dwl, 589.0e-9))
print(sig("589.6nm", 589e-9, dwl, 589.6e-9))

# Natrium
# Measured Values
wl = [331e-9, 395e-9, 404e-9, 416e-9, 420e-9, 427e-9, 430e-9, 434e-9, 439e-9, 442e-9, 449e-9, 454e-9, 467e-9, 475e-9, 498e-9, 515e-9, 525e-9, 568e-9, 589e-9, 668e-9, 675e-9, 687e-9, 696e-9, 703e-9, 707e-9, 715e-9, 720e-9, 727e-9, 738e-9, 751e-9, 763e-9, 770e-9, 773e-9, 795e-9, 801e-9, 811e-9, 819e-9]
dwl = [1.5e-9 for i in wl]

# 1. Nebenserie
ns1 = [819e-9]
m = 3
Eryd = -13.605
E3p = Eryd / m**2 - h / e * c / ns1[0]

for m in range(4,13):
  ns1.append(h / e * c / (Eryd * m**-2 - E3p))

print()
print("1. Nebenserie")
print()
print(lst("Lambda", ns1, [0.3e-9 for i in range(len(ns1))]))

# 2. Nebenserie
ns2 = 589e-9
E3s = E3p - h / e * c / ns2
sCorr = 3 - sqrt(Eryd / E3s)

ns2 = []
for m in range(4,10):
  ns2.append(h / e * c / (Eryd  * (m - sCorr)**-2 - E3p))

print()
print("2. Nebenserie")
print()
print(lst("Lambda", ns2, [0.3e-9 for i in range(len(ns2))]))

# Hauptserie
pCorr = 3 - sqrt(Eryd / E3p)

hs = []
for m in range(4,6):
  hs.append(h / e * c / (Eryd  * (m - pCorr)**-2 - E3s))

print()
print("Hauptserie")
print()
print(lst("Lambda", hs, [0.3e-9 for i in range(len(hs))]))

# Plot
initplot(title="Wellenlängen", xlabel="Wellenlänge", ylabel="Wellenlänge")
plot(ns1, ns1, label="1. Nebenserie")
plot(ns2, ns2, label="2. Nebenserie")
plot(hs, hs, label="Hauptserie")
plot(wl, wl, xerr=dwl, yerr=dwl, label="Gemessen")
#showplot()

# Identifizierungen
dwl = dwl[0]
print()
print("Abweichungen")
print()
print(sig("HS  332.1nm", wl[0], dwl, hs[0]))
print(sig("NS1 423.0nm", wl[4], dwl, ns1[9]))
print(sig("NS1 425.6nm", wl[5], dwl, ns1[8]))
print(sig("NS1 429.1nm", wl[6], dwl, ns1[7]))
print(sig("NS1 433.9nm", wl[7], dwl, ns1[6]))
print(sig("NS1 440.8nm", wl[8], dwl, ns1[5]))
print(sig("NS2 444.1nm", wl[9], dwl, ns2[5]))
print(sig("NS1 451.2nm", wl[10], dwl, ns1[4]))
print(sig("NS2 456.5nm", wl[11], dwl, ns2[4]))
print(sig("NS1 468.3nm", wl[12], dwl, ns1[3]))
print(sig("NS2 477.6nm", wl[13], dwl, ns2[3]))
print(sig("NS1 499.7nm", wl[14], dwl, ns1[2]))
print(sig("NS2 518.7nm", wl[15], dwl, ns2[2]))
print(sig("NS1 570.0nm", wl[17], dwl, ns1[1]))
print(sig("NS1 819.0nm", wl[36], dwl, ns1[0]))
