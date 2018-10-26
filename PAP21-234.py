import measure as ms

# 1. Nebenserie
wl = [819e-9]
m = 3
Eryd = -13.605
E3p = Eryd / m**2 - ms.h / ms.e * ms.c / wl[0]

for m in range(4,13):
  wl.append(ms.h / ms.e * ms.c / (Eryd * m**-2 - E3p))

print()
print("1. Nebenserie")
print()
print(ms.val("E3p", E3p))
print()
print(ms.lst("Lambda", wl, [0.3e-9 for i in range(len(wl))]))

# 2. Nebenserie
wl = 589e-9
E3s = E3p - ms.h / ms.e * ms.c / wl
sCorr = 3 - ms.sqrt(Eryd / E3s)

wl = []
for m in range(4,10):
  wl.append(ms.h / ms.e * ms.c / (Eryd  * (m - sCorr)**-2 - E3p))

print()
print("2. Nebenserie")
print()
print(ms.lst("Lambda", wl, [0.3e-9 for i in range(len(wl))]))

# Hauptserie
pCorr = 3 - ms.sqrt(Eryd / E3p)

wl = []
for m in range(4,6):
  wl.append(ms.h / ms.e * ms.c / (Eryd  * (m - pCorr)**-2 - E3s))

print()
print("Hauptserie")
print()
print(ms.lst("Lambda", wl, [0.3e-9 for i in range(len(wl))]))
