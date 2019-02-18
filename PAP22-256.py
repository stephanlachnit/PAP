# measure version 1.8.6
from measure import plt,np,npfarray,sqrt,pltext,curve_fit,val,c,h,sig

Z = npfarray([26,42,29,22,47,40,30,28])

# Literaturwerte
sqrt_Er_lit = sqrt(13.6e-3)
sig12_lit = 1
p0 = npfarray([sqrt_Er_lit,sig12_lit])

print()
print(val('Literaturwert sqrt(Er)',sqrt_Er_lit))

# K_alpha
K_alpha = npfarray([6.42,17.47,8.08,4.49,21.90,15.79,8.68,7.51])
Delta_K_alpha = npfarray([0.16,0.17,0.15,0.16,0.20,0.17,0.16,0.15])
sqrt_K_alpha = sqrt(K_alpha)
Delta_sqrt_K_alpha = 1/2 * 1/sqrt(K_alpha) * Delta_K_alpha

pltext.initplot(num=1,title=r'$\sqrt{E_\alpha}$ als Funktion von $Z$',xlabel=r'$Z$',ylabel=r'$\sqrt{E_\alpha}$ in $\sqrt{keV}$')
pltext.plotdata(x=Z, y=sqrt_K_alpha, dy=Delta_sqrt_K_alpha)

n1=1
n2=2
def fit_func_alpha(x,sqrt_Er,sig12):
  return sqrt_Er*(x-sig12)*sqrt(1/n1**2 - 1/n2**2)

popt,pcov = curve_fit(fit_func_alpha,Z,sqrt_K_alpha,sigma=Delta_sqrt_K_alpha,p0=p0)

sqrt_Er_alpha = popt[0]
Delta_sqrt_Er_alpha = sqrt(pcov[0,0])
sig12_alpha = popt[1]
Delta_sig12_alpha = sqrt(pcov[1,1])

plt.plot(Z, fit_func_alpha(Z,*popt))

print()
print('K_alpha:')
print(val('sqrt(Er)',sqrt_Er_alpha,Delta_sqrt_Er_alpha))
print(sig('Abweichung',sqrt_Er_alpha,Delta_sqrt_Er_alpha,sqrt_Er_lit,perc=True))
print(val('sig12',sig12_alpha,Delta_sig12_alpha))
print(sig('Abweichung',sig12_alpha,Delta_sig12_alpha,sig12_lit,perc=True))

# K_beta mit Ti

K_beta = npfarray([7.05,19.56,8.95,17.46,24.59,17.67,9.63,8.30])
Delta_K_beta = npfarray([0.19,0.16,0.17,0.32,0.14,0.17,0.15,0.16])
sqrt_K_beta = sqrt(K_beta)
Delta_sqrt_K_beta = 1/2 * 1/sqrt(K_beta) * Delta_K_beta

pltext.initplot(num=2,title=r'$\sqrt{E_\beta}$ als Funktion von $Z$',xlabel=r'$Z$',ylabel=r'$\sqrt{E_\beta} in \sqrt{keV}$')
pltext.plotdata(x=Z, y=sqrt_K_beta, dy=Delta_sqrt_K_beta)

n1=1
n2=3
def fit_func_beta(x,sqrt_Er,sig12):
  return sqrt_Er*(x-sig12)*sqrt(1/n1**2 - 1/n2**2)

popt,pcov = curve_fit(fit_func_beta,Z,sqrt_K_beta,sigma=Delta_sqrt_K_beta,p0=p0)

sqrt_Er_beta = popt[0]
Delta_sqrt_Er_beta = sqrt(pcov[0,0])
sig12_beta = popt[1]
Delta_sig12_beta = sqrt(pcov[1,1])

plt.plot(Z, fit_func_beta(Z,*popt))

print()
print('K_beta:')
print(val('sqrt(Er)',sqrt_Er_beta,Delta_sqrt_Er_beta))
print(sig('Abweichung',sqrt_Er_beta,Delta_sqrt_Er_beta,sqrt_Er_lit,perc=True))
print(val('sig12',sig12_beta,Delta_sig12_beta))

# K_beta ohne Ti

Z_boTi = np.delete(Z,3)
K_boTi = np.delete(K_beta,3)
Delta_K_boTi = np.delete(Delta_K_beta,3)
sqrt_K_boTi = sqrt(K_boTi)
Delta_sqrt_K_boTi = 1/2 * 1/sqrt(K_boTi) * Delta_K_boTi

pltext.initplot(num=3,title=r'$\sqrt{E_\beta}$ als Funktion von $Z$ (ohne Ti)',xlabel=r'$Z$',ylabel=r'$\sqrt{E_\beta} in \sqrt{keV}$')
pltext.plotdata(x=Z_boTi, y=sqrt_K_boTi, dy=Delta_sqrt_K_boTi)

popt,pcov = curve_fit(fit_func_beta,Z_boTi,sqrt_K_boTi,sigma=Delta_sqrt_K_boTi,p0=p0)

sqrt_Er_boTi = popt[0]
Delta_sqrt_Er_boTi = sqrt(pcov[0,0])
sig12_boTi = popt[1]
Delta_sig12_boTi = sqrt(pcov[1,1])

plt.plot(Z_boTi, fit_func_beta(Z_boTi,*popt))

print()
print('K_beta (ohne Ti):')
print(val('sqrt(Er)',sqrt_Er_boTi,Delta_sqrt_Er_boTi))
print(sig('Abweichung',sqrt_Er_boTi,Delta_sqrt_Er_boTi,sqrt_Er_lit,perc=True))
print(val('sig12',sig12_boTi,Delta_sig12_boTi))
print(sig('Abweichung',sig12_boTi,Delta_sig12_boTi,sig12_lit,perc=True))

plt.show()
