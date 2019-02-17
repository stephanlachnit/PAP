import measure as ms
from measure import pi, sqrt
import numpy as np
import scipy.constants as sc
from scipy.optimize import curve_fit

# Data of the helmholtz-coil
d_h = 0.295
h_h = 0.147
n_h = 124

# Data of the induction coil
n_i = 4000
A_i = 41.7 * sc.centi**2

# Measured values
# Resistance and maximum allowed supply voltage
I_max = 5.0
I_test = 1.83
d_I_test = 0.01
U_test = 2.2
d_U_test = 0.1
R_sp = U_test / I_test
d_R_sp = R_sp * ms.sqrt((d_U_test / U_test)**2 + (d_I_test / I_test)**2)
U_max = R_sp * I_max

# Frequency and induced voltage at constant current
I1 = 8.0 / 2
d_I1 = 0.1 / 2
f1 = [3.0, 5.6, 9.3, 11.93, 14.90]
d_f1 = [0.1, 0.1, 0.1, 0.40, 0.40]
U1_i = [0.81, 2.40, 5.08, 6.96, 9.04]
d_U1_i = [0.02, 0.03, 0.04, 0.01, 0.01]

# Current, voltage and induced voltage at constant frequency
f2 = 9.96
d_f2 = 0.05
I2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
d_I2 = [0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1]
U2 = [0.6, 1.2, 1.8, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4]
d_U2 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
U2_i = [0.855, 1.540, 2.240, 2.940, 3.560, 4.240, 4.920, 5.620, 6.260]
d_U2_i = [0.010, 0.010, 0.020, 0.020, 0.040, 0.020, 0.040, 0.020, 0.020]

# Rotation angle and induced voltage at constant a.c. voltage
f3 = 98.03
d_f3 = 1.0
alpha3 = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
d_alpha3 = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
U3_i = [2.52, 2.24, 1.36, 0.046, 1.25, 2.22, 2.51, 2.14, 1.24, 0.028, 1.19, 2.14]
d_U3_i = [0.01, 0.01, 0.01, 0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 0.002, 0.01, 0.01]

# Current, voltage, a.c. frequency and induced voltage at a.c. voltage
I4 = [0.68, 0.44, 0.31, 0.24, 0.185, 0.15, 0.13, 0.11, 0.095, 0.08, 0.0509, 0.035175, 0.02667, 0.021515, 0.01792, 0.01533, 0.013475, 0.012125, 0.01070]
d_I4 = [0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.01, 0.0001, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002, 0.00002]
U4 = [2.8, 4.20, 4.16, 4.08, 4.08, 4.04, 4.00, 4.00, 4.00, 3.96, 3.92, 3.80, 3.80, 3.80, 3.84, 3.80, 3.76, 3.80, 3.80]
d_U4 = [0.2, 0.04, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
f4 = [20.60, 40.20, 60.00, 79.87, 101.2, 123.02, 140.82, 161.13, 179.8, 200.0, 405.4, 603.26, 803.8, 1000, 1205, 1412, 1609, 1805, 2033]
d_f4 = [0.03, 0.04, 0.10, 0.20, 0.10, 0.3, 0.3, 0.3, 0.3, 0.4, 0.1, 0.10, 0.4, 2, 2, 3, 3, 3, 3]
U4_i = [3.30, 4.40, 4.76, 4.92, 5.00, 5.00, 5.04, 5.04, 5.08, 5.08, 5.08, 5.04, 5.08, 5.08, 5.16, 5.16, 5.16, 5.20, 5.32]
d_U4_i = [0.01, 0.01, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]

# Frequency and induced voltage at no current
f5 = 15.06
d_f5 = 0.10
U5_i = 0.158
d_U5_i = 0.002

# Frequency, compensating current and induced voltage
f6 = 15.1
d_f6 = 0.2
I6_c = 0.018205
d_I6_c = 0.00001
U6_i = 0.072
d_U6_i = 0.001

U1_i = [U / 2 for U in U1_i]
d_U1_i = [d_U / 2 for d_U in d_U1_i]
U2_i = [U / 2 for U in U2_i]
d_U2_i = [d_U / 2 for d_U in d_U2_i]
U3_i = [U / 2 for U in U3_i]
d_U3_i = [d_U / 2 for d_U in d_U3_i]
U4_i = [U / 2 for U in U4_i]
d_U4_i = [d_U / 2 for d_U in d_U4_i]
U5_i = U5_i / 2
d_U5_i = d_U5_i / 2
U6_i = U6_i / 2
d_U6_i = d_U6_i / 2

# Experiment 1
[s1, d_s1, tmp, tmp] = ms.linreg(f1, U1_i, d_U1_i, d_f1, plot=False, graphname='Induced voltage at constant current')
#ms.plt.show()
#ms.linreg(I2, U2_i, d_U2_i, d_I2, plot=True, graphname='Induced voltage at constant frequency')
#ms.plt.show()

B1 = s1 / (2 * pi * n_i * A_i)
d_B1 = d_s1 / (2 * pi * n_i * A_i)

B1_ = (8 / sqrt(125)) * sc.mu_0 * I1 * n_h / h_h
d_B1_ = (8 / sqrt(125)) * sc.mu_0 * d_I1 * n_h / h_h

print(ms.val("B1", B1, d_B1))
print(ms.val("B1'", B1_, d_B1_))
print(ms.sig("B1,B1'", B1, d_B1, B1_, d_B1_))
print()

# Experiment 2
#ms.pltext.initplot(title='Induktionsspannung in Abhängigkeit des Winkels', xlabel='α', ylabel='U')
#ms.pltext.plotdata(alpha3, U3_i, d_U3_i, d_alpha3, connect=True)
#ms.plt.show()

Ui_U_ratio = [U4_i[i] / U4[i] for i in range(len(U4_i))]
d_Ui_U_ratio = [Ui_U_ratio[i] * sqrt((d_U4_i[i] / U4_i[i])**2 + (d_U4[i] / U4[i])**2) for i in range(len(U4_i))]
ms.pltext.initplot(title='Verhältnis der Induktionsspannung zur Spulenspannung in Abhängigkeit der Frequenz', xlabel='f', ylabel='Ui/U')
ms.pltext.plotdata(f4, Ui_U_ratio, d_Ui_U_ratio, d_f4, connect=True)
ms.plt.show()
