import measure as ms
from measure import npfarray as npf
import numpy as np
from numpy import pi as π
from numpy import sqrt, exp
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.special import gamma

# Measured values
# Measurement of the characteristic
V_E = 440
t0_c = 30.0
U_c = npf([500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750])
n_c = npf([2462, 2606, 2636, 2738, 2753, 2671, 2618, 2685, 2742, 2715, 2792])
d_n_c = sqrt(n_c) / t0_c
n_c = n_c / t0_c

# Determination of the characteristic
#ms.pltext.initplot()
#ms.pltext.plotdata(U_c, n_c, d_n_c)
[s, d_s, i, d_i] = ms.linreg(U_c[1:], n_c[1:], d_n_c[1:])
x_lin = np.linspace(U_c[0], U_c[-1], 1000)
y_lin = s * x_lin + i
y_err_lin = (s + d_s) * x_lin + (i - d_i)
#ms.plt.plot(x_lin, y_lin, color='blue')
#ms.plt.plot(x_lin, y_err_lin, color='blue', linestyle='dashed')
#ms.plt.show()

U_0 = 600

# Measurement of the plateau-slope
t1_c = sc.minute
U1_c = npf([U_0, U_0 + 100])
n1_c = npf([12022, 12172]) / t1_c

t2_c = 3 * sc.minute
U2_c = npf([U_0, U_0 + 100])
n2_c = npf([35870, 36579]) / t2_c

# Determination of the plateau-slope
del_n1 = n1_c[1] - n1_c[0]
del_n1_p = del_n1 / n1_c[0]
d_del_n1 = sqrt(n1_c[1] - n1_c[0])
d_del_n1_p = del_n1_p * sqrt((d_del_n1 / del_n1)**2 + (sqrt(n1_c[0]) / n1_c[0])**2)

del_n2 = n2_c[1] - n2_c[0]
del_n2_p = del_n2 / n2_c[0]
d_del_n2 = sqrt(n2_c[1] - n2_c[0])
d_del_n2_p = del_n2_p * sqrt((d_del_n2 / del_n2)**2 + (sqrt(n2_c[0]) / n2_c[0])**2)

print("Plateau-slope:")
print(ms.val('del_n1', del_n1, d_del_n1))
print(ms.val('del_n1_p', del_n1_p, d_del_n1_p))
print(ms.sig('del_n1,0', del_n1, d_del_n1, 0.0))
print()
print(ms.val('del_n2', del_n2, d_del_n2))
print(ms.val('del_n2_p', del_n2_p, d_del_n2_p))
print(ms.sig('del_n2,0', del_n2, d_del_n2, 0.0))

# Statistics
def g_pdf(x, μ, σ, A):
  return A / sqrt(2 * π * σ**2) * exp(-(x - μ)**2 / (2 * σ**2))
def p_pdf(k, μ, A):
  return A * (μ**k / gamma(k + 1)) * exp(-μ)
# 1.Determination of the statistics
n3, h3 = np.loadtxt('data/PAP22/251_1j.dat', unpack=True)
d_h3 = sqrt(h3)

[μ3_g, σ3_g, A3_g], _ = curve_fit(g_pdf, n3[1:-2], h3[1:-2], sigma=d_h3[1:-2], p0=[25, 5, 2000])
[μ3_p, A3_p], _ = curve_fit(p_pdf, n3[1:-2], h3[1:-2], sigma=d_h3[1:-2], p0=[25, 2000])
x3_gp = np.linspace(n3[0], n3[-1], 1000)
y3_g = g_pdf(x3_gp, μ3_g, σ3_g, A3_g)
y3_p = p_pdf(x3_gp, μ3_p, A3_p)

ms.pltext.initplot()
ms.pltext.plotdata(n3, h3, d_h3, color='gray')
ms.plt.plot(x3_gp, y3_g)
ms.plt.plot(x3_gp, y3_p)
ms.plt.show()

# 2.Determination of the statistics
n4, h4 = np.loadtxt('data/PAP22/251_2j.dat', unpack=True)
d_h4 = sqrt(h4)

[μ4_g, σ4_g, A4_g], _ = curve_fit(g_pdf, n4[1:-2], h4[1:-2], sigma=d_h4[1:-2], p0=[4, 1.5, 4000])
[μ4_p, A4_p], _ = curve_fit(p_pdf, n4[:-2], h4[:-2], sigma=d_h4[:-2], p0=[4, 4000])
x4_gp = np.linspace(n4[0], n4[-1], 1000)
y4_g = g_pdf(x4_gp, μ4_g, σ4_g, A4_g)
y4_p = p_pdf(x4_gp, μ4_p, A4_p)

ms.pltext.initplot()
ms.pltext.plotdata(n4, h4, d_h4, color='gray')
ms.plt.plot(x4_gp, y4_g)
ms.plt.plot(x4_gp, y4_p)
ms.plt.show()
