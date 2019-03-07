# measure version 1.8.11s
from measure import npfarray,pltext,plt,sqrt,expreg

# Aufgabe 2
ug_5min = 120
ug_5min_err = sqrt(ug_5min)
ug_1s = ug_5min / (5 * 60)
ug_1s_err = ug_5min_err / (5 * 60)

# Aufgabe 3
t = npfarray([30,30,30,30,30,30,30,60,120,120,120,120,120])
d = npfarray([0,1,2,3,4,5,6,7,8,9,10,11,12]) * 0.3
N = npfarray([1196,677,449,289,189,143,78,86,123,91,61,68,65])
N_err = sqrt(N)

N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

pltext.initplot(num=1,title=r'Absoprtion $\beta$-Strahlung',xlabel='Dicke in mm',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(d,N_t_ratio,N_t_ratio_err,label='Messwerte')
pltext.set_layout()

# Aufgabe 4
d = npfarray([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0])
N = npfarray([3786,2609,1845,1376,1092,775,599,414,335,314,220])
N_err = sqrt(N)

t = 60

N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s**2)

pltext.initplot(num=2,title=r'Absorptions von $\gamma$-Strahlung (Blei)',xlabel='Dicke in cm',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(d,N_t_ratio,N_t_ratio_err,label='Messwerte')
[slope,dslope] = expreg(d,N_t_ratio,N_t_ratio_err,plot=True)[0:2]
pltext.set_layout()

# Aufgabe 5
s = npfarray([5,10,20])
s_err = npfarray([0.3,0.3,0.3])
N = npfarray([40836,11242,2889])
N_err = sqrt(N)

t = 60
N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

pltext.initplot(num=3,title=r'Absorptions von $\gamma$-Strahlung (Luft)',xlabel='Abstand in cm',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(s,N_t_ratio,N_t_ratio_err,s_err,label='Messwerte')
pltext.set_layout()

# Aufgabe 6
p = npfarray([20,57,99,126,162,222,321,353,369,391,434,470,495,521,549,602])
p_err = npfarray([2 for x in range(len(p))])

N = npfarray([4868,4803,4915,4835,4710,4689,3971,2561,1489,609,84,67,73,93,68,83])
N_err = sqrt(N)

t = 60
N_t_ratio = N / t - ug_1s
N_t_ratio_err = sqrt((N_err / t)**2 + ug_1s_err**2)

pltext.initplot(num=4,title=r'Absorption $\alpha$-Strahlung',xlabel='Druck in mbar',ylabel='Ereignisse pro Sekunde',scale='linlog')
pltext.plotdata(p,N_t_ratio,N_t_ratio_err,p_err,label='Messwerte')
pltext.set_layout(xlim=(0,625),ylim=(5e-1,1e2))

# Plot
print()
plt.show()
