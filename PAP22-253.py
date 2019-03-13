# measure version 1.9.1s
from measure import npfarray,pltext,plt,sqrt,expreg,val,linreg,tbl,lst,exp,mv,dtot_mv

# Aufgabe 2
ug_5min = 120
ug_5min_err = sqrt(ug_5min)
ug_1s = ug_5min / (5 * 60)
ug_1s_err = ug_5min_err / (5 * 60)

print('\nAufgabe 2:\n')
print(val(ug_1s,ug_1s_err,name='Untergrund / s'))

# Aufgabe 3
t = npfarray([30,30,30,30,30,30,30,60,120,120,120,120,120,120,300])
d = npfarray([0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0,3.3,3.6,3.9,4.9])
N = npfarray([1196,677,449,289,189,143,78,86,123,91,61,68,65,42,131])
N_err = sqrt(N)

N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

pltext.initplot(num=1,title=r'Abbildung   : Absoprtion $\beta$-Strahlung',xlabel='Dicke in mm',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(d,N_t_ratio,N_t_ratio_err,label='Messwerte')
pltext.set_layout(xlim=(-0.2,5.2),ylim=(1e-2,6e1))

pltext.initplot(num=2,title=r'Abbildung   : Absoprtion $\beta$-Strahlung',xlabel='Dicke in mm',ylabel='Ereignisse pro Sekunde',scale='linlin')
pltext.plotdata(d[-6:],N_t_ratio[-6:],N_t_ratio_err[-6:],label='Messwerte')
[slope,dslope,yitc,dyitc] = linreg(d[-6:],N_t_ratio[-6:],N_t_ratio_err[-6:],plot=True,prange=(2.5,5))
pltext.set_layout(xlim=(2.5,5),ylim=(-0.15,0.5))

from scipy.optimize import fsolve
def lin(x):
  return slope * x + yitc
def lin_err(x):
  return (slope - dslope) * x + (yitc + dyitc)

dmax = fsolve(lin,x0=4.8)[0] * 1e-3
dmax_err = dmax - fsolve(lin_err,x0=4.2)[0] * 1e-3

dichte = 2.70e3
dichte_err = 0.01e3
R_ES = 0.13e1

Rbeta = dichte * dmax + R_ES 
Rbeta_err = sqrt((dichte * dmax_err)**2 + (dichte_err * dmax)**2)

E = 2.7e6

print('\nAufgabe 3:\n')
print(val(dmax,dmax_err,name='dmax'))
print(val(Rbeta,Rbeta_err,name='Rbeta'))
print(val(E,name='E'))

# Aufgabe 4
d = npfarray([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])*1e1
N = npfarray([3786,2609,1845,1376,1092,775,599,414,335,314,220])
N_err = sqrt(N)

t = 60

N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s**2)

pltext.initplot(num=3,title=r'Abbildung   : Absorptions von $\gamma$-Strahlung (Blei)',xlabel='Dicke in mm',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(d,N_t_ratio,N_t_ratio_err,label='Messwerte')
[slope,slope_err] = expreg(d,N_t_ratio,N_t_ratio_err,plot=True,prange=(-2.5,52.5))[0:2]
pltext.set_layout(xlim=(-2.5,52.5),ylim=(2e0,7e1))

mu = -slope * 1e3
mu_err = slope_err * 1e3
roh = 11.34e3
roh_err = 0.01e3
msk = mu / roh
msk_err = 1 / roh * sqrt(mu_err**2 + (mu / roh * roh_err)**2)
E = 1.3e6

print('\nAufgabe 4:\n')
print(val(mu,mu_err,name='mu'))
print(val(msk,msk_err,name='msk'))
print(val(E,name='E'))

# Aufgabe 5
s = npfarray([5,10,20])*1e-2
s_err = npfarray([0.3,0.3,0.3])*1e-2
N = npfarray([40836,11242,2889])
N_err = sqrt(N)

t = 60
N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

radius = 7e-3
epsilon = 0.04
l = 4e-2

A = 4 * N_t_ratio * s**2 / (epsilon * radius**2)
A_err =  4 * s / (epsilon * radius**2) * sqrt((s * N_t_ratio_err)**2 + (N_t_ratio * 2 * s_err)**2)

k1 = (s + l/2)**2 / s**2
k1_err = l/s**2 * (1 + l/(2 * s)) * s_err

A1 = k1 * A
A1_err = sqrt((k1_err * A)**2 + (k1 * A_err)**2)

d = 1.4e-3
dichte = 7.9e3
k2 = exp(-msk * dichte * d)
k2_err = dichte * d * msk_err * exp(-msk * dichte * d)

A2 = k2 * A1
A2_err = sqrt((k2_err * A1)**2 + (k2 * A1_err)**2)

print('\nAufgabe 5:\n')
print(tbl([lst(s,s_err,name='s / m'),lst(A,A_err,name='A / Bq'),lst(A1,A1_err,name='A1 / Bq'),lst(A2,A2_err,name='A2 / Bq')]))

# Aufgabe 6
p = npfarray([20,57,99,126,162,222,321,353,369,391,434,470,495,521,549,602])
p_err = npfarray([2 for x in range(len(p))])

N = npfarray([4868,4803,4915,4835,4710,4689,3971,2561,1489,609,84,67,73,93,68,83])
N_err = sqrt(N)

t = 60
N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

pltext.initplot(num=4,title=r'Abbildung   : Absorption $\alpha$-Strahlung',xlabel='Druck in mbar',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(p,N_t_ratio,N_t_ratio_err,p_err,label='Messwerte')
expreg(p[8:11],N_t_ratio[8:11],N_t_ratio_err[8:11],plot=True,prange=(0,625))
pltext.set_layout(xlim=(0,625),ylim=(5e-1,1e2))

p_halb = 391
p_halb_err = 8
p0 = 1013
s0 = 3.95e-2
s0_err = 0.05e-2
roh_glf = 2.35
roh_brems = 1.43

s1 = mv(p / p0 * s0)
s1_err = dtot_mv(p / p0 * s0,1/p0 * sqrt((p_err * s0)**2 + (p * s0_err)**2))

s2 = s1 + roh_glf / roh_brems * 1e-2
s2_err = s1_err

s3 = s2 + 0.68e-2
s3_err = s2_err

E = 5.3e6

print('\nAufgabe 6:\n')
print(val(p_halb,p_halb_err,name='p_halb'))
print()
print(tbl([lst([s1],[s1_err],name='s1 / m'),lst([s2],[s2_err],name='s2 / m'),lst([s3],[s3_err],name='s3 / m')]))
print()
print(val(E,name='E'))

# Plot
print()
plt.show()
