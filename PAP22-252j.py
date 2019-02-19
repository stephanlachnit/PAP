import measure as ms
from measure import sqrt, exp
import numpy as np
import scipy.constants as cs
from scipy.optimize import curve_fit

# Measured data
U0 = 520
d_U0 = 10
T_tor = 10.0
T0 = 8 * cs.minute
T_Ag = 400
T_In = 50 * cs.minute

n0 = np.loadtxt('data/PAP22/252_1j.dat', usecols=[1], unpack=True)
n_Ag_1 = np.loadtxt('data/PAP22/252_2j.dat', usecols=[1], unpack=True)
n_Ag_2 = np.loadtxt('data/PAP22/252_3j.dat', usecols=[1], unpack=True)
n_Ag_3 = np.loadtxt('data/PAP22/252_4j.dat', usecols=[1], unpack=True)
n_Ag_4 = np.loadtxt('data/PAP22/252_5j.dat', usecols=[1], unpack=True)
n_In = np.loadtxt('data/PAP22/252_6j.dat', usecols=[1], unpack=True)

# Background radiation
n0_m = ms.mv(n0)
d_n0_m = ms.dsto_mv(n0)

# Ag decay
t_Ag = np.arange(T_tor / 2, T_Ag + T_tor / 2, T_tor)
N_Ag = n_Ag_1 + n_Ag_2 + n_Ag_3 + n_Ag_4
d_N_Ag = sqrt(N_Ag)

def f(x, A1, λ1, A2, λ2):
  return A1 * exp(-λ1 * x) + A2 * exp(-λ2 * x) + 4 * n0_m

[A1, λ1, A2, λ2], pcov = curve_fit(f, t_Ag, N_Ag, sigma=d_N_Ag, p0=[250, 0.03, 250, 0.03])
[d_A1, d_λ1, d_A2, d_λ2] = np.diag(pcov)

print(ms.val("A1", A1, d_A1))
print(ms.val("λ1", λ1, d_λ1))
print(ms.val("A2", A2, d_A2))
print(ms.val("λ2", λ2, d_λ2))

t1 = np.linspace(t_Ag[0], t_Ag[-1], 1000)
N1 = f(t1, A1, λ1, A2, λ2)
ms.pltext.initplot(num=1)
ms.pltext.plotdata(t_Ag, N_Ag, d_N_Ag, color='gray')
ms.plt.plot(t1, N1)

# In decay

ms.plt.show()
