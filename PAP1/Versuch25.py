import math as m

T = [8.07e-6, 5.64e-6, 43.4e-6, 13.2e-6, 43.1e-6]
d_T = [0.03e-6, 0.02e-6, 0.02e-6, 0.1e-6, 0.1e-6]
Uss = [2.36, 10.5, 0.584, 3.54, 24.2]
d_Uss = [0.02, 0.2, 0.005, 0.02, 0.02]
Udc = [10e-3, 0.0, 2e-3, 10e-3, 100e-3]
d_Udc = [0.0, 0.0, 0.0, 0.0, 0.0]
Umax = [0.0, 0.0, 0.0, 0.0, 0.0]
d_Umax = [0.0, 0.0, 0.0, 0.0, 0.0]
Umin = [0.0, 0.0, 0.0, 0.0, 0.0]
d_Umin = [0.0, 0.0, 0.0, 0.0, 0.0]
for i in range(5):
 d_T[i] = m.sqrt(d_T[i]**2 + (0.03 * T[i])**2)
 d_Uss[i] = m.sqrt(d_Uss[i]**2 + (0.03 * Uss[i])**2)
 d_Udc[i] = 0.03 * Udc[i]
 Umax[i] = Uss[i] / 2 + Udc[i]
 d_Umax[i] = m.sqrt(d_Uss[i]**2 / 4 + d_Udc[i]**2) 
 Umin[i] = - Uss[i] / 2 + Udc[i]
 d_Umin[i] = d_Umax[i]

f = [0.0, 0.0, 0.0, 0.0, 0.0]
d_f = [0.0, 0.0, 0.0, 0.0, 0.0]
for i in range(5):
 f[i] = 1/T[i]
 d_f[i] = d_T[i] / T[i]**2 

thalb = 3.0e-6
d_thalb = m.sqrt((0.17e-6)**2 + (0.03 * thalb)**2)

f = 10e3
d_f = 0.03 * f
t = [35.5e-6, 8.9e-6]
d_t = [0.5e-6, 0.5e-6]
phiyt = [0.0, 0.0]
d_phiyt = [0.0, 0.0]
b = [401e-3, 260e-3]
d_b = [2e-3, 2e-3]
a = [502e-3, 500e-3]
d_a = [4e-3, 4e-3]
phixy = [0.0, 0.0]
d_phixy = [0.0, 0.0]
for i in range(2):
 d_t[i] = m.sqrt(d_t[i]**2 + (0.03 * t[i])**2)
 d_a[i] = m.sqrt(d_a[i]**2 + (0.03 * a[i])**2)
 d_b[i] = m.sqrt(d_b[i]**2 + (0.03 * b[i])**2)
 phiyt[i] = f * t[i] * 360.0
 d_phiyt[i] = 360.0 * m.sqrt(f**2 * d_t[i]**2 + t[i]**2 * d_f**2)
 phixy[i] = m.asin(b[i] / a[i]) * 180.0 / m.pi
 d_phixy[i] = 1.0 / (a[i]**2 * abs(m.cos(phixy[i]))) * m.sqrt(a[i]**2 * d_b[i]**2 + b[i]**2 * d_a[i]**2) * 180 / m.pi
 phixy[0] = 180.0 - phixy[0]

#print()
#print(T[4]*1e6)
#print(d_T[4]*1e6)
#print()
#print(f[4]*1e-3)
#print(d_f[4]*1e-3)
#print()
#print(Uss[4])
#print(d_Uss[4])
#print()
#print(Udc[4]*1e3)
#print(d_Udc[4]*1e3)
#print()
#print(Umax[4])
#print(d_Umax[4])
#print()
#print(Umin[4])
#print(d_Umin[4])
#print()
#print(thalb*1e6)
#print(d_thalb*1e6)
print()
print(phiyt)
print(d_phiyt)
print()
print(phixy)
print(d_phixy)
