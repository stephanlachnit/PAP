import math as m
import measure as ms

def mtd(deg, mins):
  return deg + mins / 60.0

lamHg = [690.7, 623.4, 579.1, 577.0, 546.1, 499.2, 491.6, 435.8, 407.8, 404.7]

phiHg = [mtd(181, 56), mtd(180, 34), mtd(180, 18), mtd(180, 16), mtd(179, 59), mtd(179, 30), mtd(179, 26), mtd(178, 31), mtd(177, 55), mtd(177, 52)]
dphiHg = [mtd(0, 5)] * len(phiHg)
phiHe = [mtd(180, 49), mtd(180, 22), mtd(179, 34), mtd(179, 28), mtd(179, 9), mtd(178, 48)]
dphiHe = [mtd(0, 5)] * len(phiHe)
phiH = [mtd(180, 45), mtd(179, 21), mtd(178, 29), mtd(178, 0)]
dphiH = [mtd(0, 5)] * len(phiH)

for i in range(len(phiHg)):
  ms.pv("lambda[" + str(i) + "]", lamHg[i])
print()
for i in range(len(phiHg)):
  ms.pve("deltaHg(lambda)[" + str(i) + "]", phiHg[i], dphiHg[i])
print()
for i in range(len(phiHe)):
  ms.pve("deltaHe(lambda)[" + str(i) + "]", phiHe[i], dphiHe[i])
print()
for i in range(len(phiH)):
  ms.pve("deltaH(lambda)[" + str(i) + "]", phiH[i], dphiH[i])


for i in range(len(phiHg)):
  phiHg[i] *= m.pi / 180.0
for i in range(len(phiHe)):
  phiHe[i] *= m.pi / 180.0
for i in range(len(phiH)):
  phiH[i] *= m.pi / 180.0
 
