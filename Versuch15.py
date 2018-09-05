import math as m
import measure as ms

# Erbeschleunigung
g = 9.80984
d_g = 2.0e-5
# Vollzylinder
r = 25.0e-3
d_r = 0.025e-3
m_vz = 445.1e-3
d_m_vz = 0.1e-3
# Hohlzylinder
r2 = 25.0e-3
d_r2 = 0.025e-3
r12 = 4.6e-3
d_r12 = 0.1e-3
r1 = r2 - r12
d_r1 = m.sqrt(d_r2**2 + d_r12**2)
m_hz = 445.2e-3
d_m_hz = 0.1e-3
# schiefe Ebene
h2 = 105.5e-3
d_h2 = 0.5e-3
h1 = 18.0e-3
d_h1 = 0.5e-3
s = 872.0e-3
d_s = 1.0e-3
alpha = m.asin((h2 - h1) / s)
d_alpha = m.sqrt( (d_h1**2 + d_h2**2 + ((h2 - h1) * d_s / s)**2) / (s**2 - (h2 - h1)**2))
h = r * (m.sin(alpha) + m.cos(alpha) - 1) + h2
d_h = m.sqrt(((m.sin(alpha) + m.cos(alpha) - 1) * d_r)**2 + d_h2**2 + (r * (m.cos(alpha) - m.sin(alpha)) * d_alpha)**2)
# Lichtschranken
s1 = 40.0e-3
d_s1 = 3.0e-3
s2 = 160.0e-3
d_s2 = 3.0e-3
s3 = 360.0e-3
d_s3 = 3.0e-3
s4 = 640.0e-3
d_s4 = 3.0e-3
s5 = 60.0e-3
d_s5 = 1.0e-3
s6 = 101e-3
d_s6 = 1.0e-3
s56 = s6 - s5
d_s56 = m.sqrt(d_s6**2 + d_s5**2)
# Zeiten
tmp = [0.360, 0.365, 0.369, 0.357, 0.361]
t1_vz = ms.mean_value(tmp)
d_t1_vz = ms.std_dev_m(tmp)
tmp = [0.725, 0.730, 0.734, 0.722, 0.725]
t2_vz = ms.mean_value(tmp)
d_t2_vz = ms.std_dev_m(tmp)
tmp = [1.085, 1.091, 1.093, 1.081, 1.085]
t3_vz = ms.mean_value(tmp)
d_t3_vz = ms.std_dev_m(tmp)
tmp = [1.440, 1.445, 1.447, 1.435, 1.438]
t4_vz = ms.mean_value(tmp)
d_t4_vz = ms.std_dev_m(tmp)
tmp = [1.744, 1.746, 1.741, 1.744, 1.740]
t5_vz = ms.mean_value(tmp)
d_t5_vz = ms.std_dev_m(tmp)
tmp = [1.780, 1.782, 1.777, 1.781, 1.777]
t6_vz = ms.mean_value(tmp)
d_t6_vz = ms.std_dev_m(tmp)
t56_vz = t6_vz - t5_vz
d_t56_vz = m.sqrt(d_t6_vz**2 + d_t5_vz**2)
tmp = [0.382, 0.407, 0.390, 0.390, 0.386]
t1_hz = ms.mean_value(tmp)
d_t1_hz = ms.std_dev_m(tmp)
tmp = [0.787, 0.813, 0.796, 0.791, 0.788]
t2_hz = ms.mean_value(tmp)
d_t2_hz = ms.std_dev_m(tmp)
tmp = [1.186, 1.214, 1.197, 1.189, 1.187]
t3_hz = ms.mean_value(tmp)
d_t3_hz = ms.std_dev_m(tmp)
tmp = [1.580, 1.607, 1.588, 1.579, 1.579]
t4_hz = ms.mean_value(tmp)
d_t4_hz = ms.std_dev_m(tmp)
tmp = [1.939, 1.926, 1.929, 1.931, 1.934]
t5_hz = ms.mean_value(tmp)
d_t5_hz = ms.std_dev_m(tmp)
tmp = [1.976, 1.967, 1.965, 1.972, 1.975]
t6_hz = ms.mean_value(tmp)
d_t6_hz = ms.std_dev_m(tmp)
t56_hz = t6_hz - t5_hz
d_t56_hz = m.sqrt(d_t6_hz**2 + d_t5_hz**2)
# Beschleunigung
a_vz = 2 / 3 * g * m.sin(alpha)
d_a_vz = 2 / 3 * m.sqrt((m.sin(alpha) * d_g)**2 + (g * m.cos(alpha) * d_alpha)**2)
a_hz = (2 * r2**2 * g * m.sin(alpha)) / (r1**2 + 3 * r2**2)
d_a_hz = (2 * r2) / (r1**2 + 3 * r2**2) * m.sqrt(r2**2 * ((m.sin(alpha) * d_g)**2 + (g * m.cos(alpha) * d_alpha)**2) + ((2 * r1 * g * m.sin(alpha)) / (r1**2 + 3 * r2**2))**2 * ((r2 * d_r1)**2 + (r1 * d_r2)**2))
# Geschwindigkeit
v_vz = s56 / t56_vz
d_v_vz = m.sqrt(d_s56**2 + (s56 * d_t56_vz / t56_vz)**2) / t56_vz
v_hz = s56 / t56_hz
d_v_hz = m.sqrt(d_s56**2 + (s56 * d_t56_hz / t56_hz)**2) / t56_hz
# Energie
e1_vz = m_vz * g * h
d_e1_vz = m.sqrt((g * h * d_m_vz)**2 + (m_vz * h * d_g)**2 + (m_vz * g * d_h)**2)
e2_vz = 3 / 4 * m_vz * v_vz**2
d_e2_vz = 3 / 4 * v_vz * m.sqrt((v_vz * d_m_vz)**2 + 4 * (m_vz * d_v_vz)**2)
e1_hz = m_hz * g * h
d_e1_hz = m.sqrt((g * h * d_m_hz)**2 + (m_hz * h * d_g)**2 + (m_hz * g * d_h)**2)
e2_hz = m_hz * v_hz**2 / 4 * (3 + (r1 / r2)**2)
d_e2_hz = v_hz / 2 * m.sqrt((3 + (r1 / r2))**2 * ((v_hz * d_m_hz)**2 + (m_hz * d_v_hz)**2) + (m_hz * v_hz * r1 / r2)**2 * (d_r1**2 + (r1 * d_r2 / r2)**2))

#print()
#print(alpha*180/m.pi)
#print(d_alpha*180/m.pi)
print()
print(h)
print(d_h)
#print()
#print(t1_vz)
#print(d_t1_vz)
#print(t2_vz)
#print(d_t2_vz)
#print(t3_vz)
#print(d_t3_vz)
#print(t4_vz)
#print(d_t4_vz)
#print()
#print(t5_vz)
#print(d_t5_vz)
#print(t6_vz)
#print(d_t6_vz)
#print(t56_vz)
#print(d_t56_vz)
#print()
#print(v_vz)
#print(d_v_vz)
#print()
#print(a_vz)
#print(d_a_vz)
print()
print(e1_vz)
print(d_e1_vz)
print(e2_vz)
print(d_e2_vz)
#print()
#print(t1_hz)
#print(d_t1_hz)
#print(t2_hz)
#print(d_t2_hz)
#print(t3_hz)
#print(d_t3_hz)
#print(t4_hz)
#print(d_t4_hz)
#print()
#print(t5_hz)
#print(d_t5_hz)
#print(t6_hz)
#print(d_t6_hz)
#print(t56_hz)
#print(d_t56_hz)
#print()
#print(v_hz)
#print(d_v_hz)
#print()
#print(a_hz)
#print(d_a_hz)
print()
print(e1_hz)
print(d_e1_hz)
print(e2_hz)
print(d_e2_hz)
