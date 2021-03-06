# measure version 1.8.12s
from measure import plt,pltext,npfarray,linreg,sqrt,val,dev,exp,pi,spcurvefit,nplinspace,tbl
from numpy import loadtxt
from scipy.special import gamma

# Aufgabe 2
U = npfarray([420,445,470,495,520,545,570])
N = npfarray([1887,2330,2337,2359,2407,2374,2310])
N_dsto = sqrt(N)

pltext.initplot(num=1,title='Abbildung   : Zählrohrcharakteristik',xlabel='Spannung in V',ylabel='Ereignisse')
pltext.plotdata(U,N,N_dsto,label='Messwerte')
[slope,dslope] = linreg(U[1:],N[1:],N_dsto[1:],plot=True,prange=(410,580))[0:2]
pltext.set_layout(xlim=(410,580),ylim=(1800,2500))

U0 = 510

print('\nAufgabe 2\n')
print(val(slope,dslope,name='Steigung'))

# Aufgabe 3
U = [510,610]
N1min = npfarray([9838,9871])
N1min_dsto = sqrt(N1min)
N3min = npfarray([29505,30144])
N3min_dsto = sqrt(N3min)

anstieg_1min = N1min[1] - N1min[0]
anstieg_1min_dsto = sqrt(N1min_dsto[1]**2 + N1min_dsto[0]**2)
rel_anstieg_1min = anstieg_1min / N1min[0]
rel_anstieg_1min_dsto = 1/N1min[0] * sqrt(anstieg_1min_dsto**2 + (anstieg_1min * N1min_dsto[0] / N1min[0])**2)

anstieg_3min = N3min[1] - N3min[0]
anstieg_3min_dsto = sqrt(N3min_dsto[1]**2 + N3min_dsto[0]**2)
rel_anstieg_3min = anstieg_3min / N3min[0]
rel_anstieg_3min_dsto = 1/N3min[0] * sqrt(anstieg_3min_dsto**2 + (anstieg_3min * N3min_dsto[0] / N3min[0])**2)

def fitfunc_3(t,c):
  return c*t

p_opt_0,p_err_0 = spcurvefit(fitfunc_3,[60,180],[N1min[0],N3min[0]],[N1min_dsto[0],N3min_dsto[0]])
p_opt_100,p_err_100 = spcurvefit(fitfunc_3,[60,180],[N1min[1],N3min[1]],[N1min_dsto[1],N3min_dsto[1]])
k0 = p_opt_0[0]
k0_err = p_err_0[0]
k100 = p_opt_100[0]
k100_err = p_err_100[0]
t_1perc = 10**4 * k100 / k0**2 * (1 + k100 / k0)
k_delta = k100 - k0
k_delta_err = sqrt(k100_err**2 + k0_err**2)

t_array = nplinspace(50,190)
pltext.initplot(num=2,title='Abbildung   : Ereignisse als Funktion der Zeit',xlabel='Zeit in s',ylabel='Ereignisse')
pltext.plotdata([60,180],[N1min[0],N3min[0]],[N1min_dsto[0],N3min_dsto[0]],label=r'Messwerte $U_0+0V$')
pltext.plotdata([60,180],[N1min[1],N3min[1]],[N1min_dsto[1],N3min_dsto[1]],label=r'Messwerte $U_0+100V$')
plt.plot(t_array,fitfunc_3(t_array,*p_opt_0),label=r'Fit $U_0 +0V$')
plt.plot(t_array,fitfunc_3(t_array,*p_opt_100),label=r'Fit $U_0 +100V$')
pltext.set_layout(xlim=(50,190),ylim=(0.9e4,3.1e4))

print('\nAufgabe 3\n')
print(val(anstieg_1min,anstieg_1min_dsto,name='abs Anstieg (1min)',percerr=True))
print(val(rel_anstieg_1min,rel_anstieg_1min_dsto,name='rel Anstieg (1min)',percerr=True))
print(val(anstieg_3min,anstieg_3min_dsto,name='abs Anstieg (3min)',percerr=True))
print(val(rel_anstieg_3min,rel_anstieg_3min_dsto,name='rel Anstieg (3min)',percerr=True))
print(dev(anstieg_1min,anstieg_1min_dsto,anstieg_3min,anstieg_3min_dsto,'Abw abs Anstiege'))
print(dev(rel_anstieg_1min,rel_anstieg_1min_dsto,rel_anstieg_3min,rel_anstieg_3min_dsto,'Abw rel Anstiege'))
print(val(k0,k0_err,name='Zählrate U0+0V',percerr=True))
print(val(k100,k100_err,name='Zählrate U0+100V',percerr=True))
print(val(k_delta,k_delta_err,name='Zählrate U+100V',percerr=True))
print(val(t_1perc,name='t (1%)'))

# Aufgabe 4
anzahl,ereig = loadtxt('data/251_A4.dat',unpack=True)
fehler = sqrt(ereig)

def gauss(x,A,mu,sig):
  return A / (sqrt(2*pi) * sig) * exp(-(x-mu)**2  / (2 * sig**2))

def poisson(x,A,mu):
  return A * exp(-mu) * mu**x / gamma(x+1)

p_opt_g,p_err_g = spcurvefit(gauss,anzahl[5:-12],ereig[5:-12],yerr=fehler[5:-12],p0=[2300,70,8])
p_opt_p,p_err_p = spcurvefit(poisson,anzahl[5:-12],ereig[5:-12],yerr=fehler[5:-12],p0=[2300,70])

sigma_theo = sqrt(p_opt_g[1])
sigma_theo_err = 1 / (2 * sqrt(p_opt_g[1])) * p_err_g[1]

from numpy import sum
from scipy.stats import chi2
chi2_g = sum((gauss(anzahl[5:-12],*p_opt_g) - ereig[5:-12])**2 / fehler[5:-12]**2)
dof_g = len(anzahl[5:-12]) - 3
chi2_red_g = chi2_g / dof_g
chi2_p = sum((poisson(anzahl[5:-12],*p_opt_p) - ereig[5:-12])**2 / fehler[5:-12]**2)
dof_p = len(anzahl[5:-12]) - 2
chi2_red_p = chi2_p / dof_p
prob_g = 1 - chi2.cdf(chi2_g,dof_g)
prob_p = 1 - chi2.cdf(chi2_p,dof_p)

x_array = nplinspace(40,100)
pltext.initplot(num=3,title='Abbildung   : Statistik des radioaktiven Zerfalls',xlabel='Zerfälle in 1/s',ylabel='Häufigkeit')
pltext.plotdata(anzahl,ereig,fehler,label='Messwerte')
plt.plot(x_array,gauss(x_array,*p_opt_g),label='Gauss Fit')
plt.plot(x_array,poisson(x_array,*p_opt_p),label='Poisson Fit')
pltext.set_layout(xlim=(40,100),ylim=(-2,142))

print('\nAufgabe 4:\n')
print(tbl([['Fitmethode','A','mu','sig'],['Gauss',val(p_opt_g[0],p_err_g[0]),val(p_opt_g[1],p_err_g[1]),val(p_opt_g[2],p_err_g[2])],['Poisson',val(p_opt_p[0],p_err_p[0]),val(p_opt_p[1],p_err_p[1]),'']]))
print()
print(val(sigma_theo,sigma_theo_err,name='sig (theo)'))
print(dev(p_opt_g[2],p_err_g[2],sigma_theo,sigma_theo_err,name='Abweichung'))
print()
print(tbl([['Fitmethode','Chi2','Chi2red','Fitwkeit'],['Gauss',val(chi2_g),val(chi2_red_g),val(prob_g)],['Poisson',val(chi2_p),val(chi2_red_p),val(prob_p)]]))

# Aufgabe 5
anzahl,ereig = loadtxt('data/251_A5.dat',unpack=True)
fehler = sqrt(ereig)

p_opt_g,p_err_g = spcurvefit(gauss,anzahl[0:-2],ereig[0:-2],yerr=fehler[0:-2],p0=[5200,4,2])
p_opt_p,p_err_p = spcurvefit(poisson,anzahl[0:-2],ereig[0:-2],yerr=fehler[0:-2],p0=[5200,4])

sigma_theo = sqrt(p_opt_g[1])
sigma_theo_err = 1 / (2 * sqrt(p_opt_g[1])) * p_err_g[1]

chi2_g = sum((gauss(anzahl[0:-2],*p_opt_g) - ereig[0:-2])**2 / fehler[0:-2]**2)
dof_g = len(anzahl[0:-2]) - 3
chi2_red_g = chi2_g / dof_g
chi2_p = sum((poisson(anzahl[0:-2],*p_opt_p) - ereig[0:-2])**2 / fehler[0:-2]**2)
dof_p = len(anzahl[0:-2]) - 2
chi2_red_p = chi2_p / dof_p
prob_g = 1 - chi2.cdf(chi2_g,dof_g)
prob_p = 1 - chi2.cdf(chi2_p,dof_p)

x_array = nplinspace(0,14)
pltext.initplot(num=4,title='Abbildung   : Statistik des radioaktiven Zerfalls',xlabel='Zerfälle in 1/s',ylabel='Häufigkeit',scale='linlog')
pltext.plotdata(anzahl,ereig,fehler,label='Messwerte')
plt.plot(x_array,gauss(x_array,*p_opt_g),label='Gauss Fit')
plt.plot(x_array,poisson(x_array,*p_opt_p),label='Poisson Fit')
pltext.set_layout(xlim=(0,14),ylim=(4e-1,2e3))

print('\nAufgabe 5:\n')
print(tbl([['Fitmethode','A','mu','sig'],['Gauss',val(p_opt_g[0],p_err_g[0]),val(p_opt_g[1],p_err_g[1]),val(p_opt_g[2],p_err_g[2])],['Poisson',val(p_opt_p[0],p_err_p[0]),val(p_opt_p[1],p_err_p[1]),'']]))
print()
print(val(sigma_theo,sigma_theo_err,name='sig (theo)'))
print(dev(p_opt_g[2],p_err_g[2],sigma_theo,sigma_theo_err,name='Abweichung'))
print()
print(tbl([['Fitmethode','Chi2','Chi2red','Fitwkeit'],['Gauss',val(chi2_g),val(chi2_red_g),val(prob_g)],['Poisson',val(chi2_p),val(chi2_red_p),val(prob_p)]]))

# Plot
print()
plt.show()
