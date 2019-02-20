# measure version 1.8.4
from measure import plt,pltext,np,npfarray,linreg,sqrt,val,curve_fit

laser = 635e-9
schirmabstand = 1

# Eichung
links = npfarray([554.48,584.99,616.35,644.31,677.36])
rechts = npfarray([446.01,413.80,383.30,355.33,322.28])
diff = links - rechts
diff_err = sqrt(2)*3 *diff/diff
sbreite = 2e-3/100*npfarray([43,70,90,104,132])
sbreite_err = 1e-3/100 *sbreite/sbreite
pltext.initplot(num=0, title='Eichung',xlabel='Pixel',ylabel='Meter')
[mpp, dmpp, itc, ditc] = linreg(x=diff,y=sbreite,dx=diff_err,dy=sbreite_err,plot=True)
print()
print(val('Meter pro Pixel',mpp,dmpp))

# Einzelspaltbreite
lage_es_hmax = 499.79
lage_es_min = (npfarray([529.91,561.26,590.92,623.13,653.63]) - lage_es_hmax)*mpp
lage_es_max = (lage_es_hmax - npfarray([456.82,424.68,394.69,363.78,333.80]))*mpp
lage_es_d = npfarray([3,3,3,3,3])*mpp
n_es = npfarray([1,2,3,4,5])

pltext.initplot(num=1,title='Position Extrema am Einzelspalt',xlabel='Ordnung',ylabel='Abstand')
[es_sl, es_dsl, es_itc, es_ditc] = linreg(n_es,lage_es_min,dy=lage_es_d)
plt.plot([0,6],[es_itc,es_itc+6*es_sl])
pltext.plotdata(n_es,lage_es_min,dy=lage_es_d,label='Minima')
pltext.plotdata((lage_es_max-es_itc)/es_sl,lage_es_max,dx=lage_es_max/es_sl+n_es*(es_dsl/es_sl-1)-(es_itc+es_ditc)/es_sl,label='Maxima')
plt.xlim((0.8,6))
plt.ylim((0.0,0.0014))
plt.legend(loc=2)

es_breite = schirmabstand * laser / es_sl
es_breite_err = schirmabstand * laser * es_dsl / es_sl**2
print()
print(val('ES Breite',es_breite,es_breite_err))

# Doppelspaltabstand
lage_ds_hmax = 504.81
lage_ds_max = (lage_ds_hmax - npfarray([492.78,483.15,474.73,464.80,454.88]))*mpp
lage_ds_min = (lage_ds_hmax - npfarray([498.20,485.26,477.14,469.62,459.39]))*mpp
lage_ds_d = npfarray([3,3,3,3,3])*mpp
n_ds = npfarray([1,2,3,4,5])

pltext.initplot(num=2,title='Position Extrema am Doppelspalt',xlabel='Ordnung',ylabel='x in m')
[ds_sl, ds_dsl, ds_itc, ds_ditc] = linreg(n_ds,lage_ds_min,dy=lage_ds_d)
plt.plot([0,6],[ds_itc,ds_itc+6*ds_sl])
pltext.plotdata(n_ds,lage_ds_min, dy=lage_ds_d,label='Minima')
pltext.plotdata((lage_ds_max-ds_itc)/ds_sl,lage_ds_max,dx=lage_ds_max/ds_sl+n_ds*(ds_dsl/ds_sl-1)-(ds_itc+ds_ditc)/ds_sl,label='Maxima')
plt.xlim((0.8,6))
plt.ylim((0.0,0.0004))
plt.legend(loc=2)

ds_abstand = schirmabstand * laser / ds_sl
ds_abstand_err = schirmabstand * laser * ds_dsl / ds_sl**2
print()
print(val('DS Abstand',ds_abstand,ds_abstand_err))

# Intensität ES
I_es_hmax = 3840
I_es_max = [299,139,139/3827 * 2355, 139/3827 * 1586, 139/3827 * 1018]
I_es = npfarray(I_es_max[::-1] + [I_es_hmax] + I_es_max)/I_es_hmax
I_es_err = npfarray([139/3827*5,139/3827*5,139/3827*5,5,5,5,5,5,139/3827*5,139/3827*5,139/3827*5])/I_es_hmax
x_es = np.concatenate((-lage_es_max[::-1],[0],lage_es_max))

x_es_min = np.concatenate((-lage_es_min[::-1],lage_es_min))
I_es_min = 0* x_es_min

pltext.initplot(num=3, title='Intensitätsverteilung ES', xlabel='x in m', ylabel='Intensität')
pltext.plotdata(x_es,I_es,I_es_err,label='Maxima')
pltext.plotdata(x_es_min,I_es_min,label='Minima')

def I_es_theo(x):
  alpha = np.arctan2(x,schirmabstand)
  x_skr = np.pi * es_breite * np.sin(alpha) / laser
  return (np.sin(x_skr)/x_skr)**2

x_es_theo = np.linspace(-1.3e-3,1.3e-3,2000)
plt.plot(x_es_theo,I_es_theo(x_es_theo),label='berechnet')
plt.legend()

# Intensität DS
I_ds_hmax = 4032
I_ds_max = [3005,3005/3918*904,3005/3918*450,3005/3918*359,3005/3918*293]
I_ds = npfarray(I_ds_max[::-1] + [I_ds_hmax] + I_ds_max)/I_ds_hmax
I_ds_err = npfarray([3005/3918*5,3005/3918*5,3005/3918*5,3005/3918*5,5,5,5,3005/3918*5,3005/3918*5,3005/3918*5,3005/3918*5])/I_ds_hmax
x_ds = np.concatenate((-lage_ds_max[::-1],[0],lage_ds_max))

x_ds_min = np.concatenate((-lage_ds_min[::-1],lage_ds_min))
I_ds_min = 0* x_ds_min

pltext.initplot(num=4, title='Intensitätsverteilung DS', xlabel='x in m', ylabel='Intensität')
pltext.plotdata(x_ds,I_ds,I_ds_err,label='Maxima')
pltext.plotdata(x_ds_min,I_ds_min,label='Minima')

max1 = [lage_ds_max[0],I_ds_max[0]/I_ds_hmax]
from scipy.optimize import fsolve
def ds_breite_nst(b):
  alpha = np.arctan2(max1[0],schirmabstand)
  x_skr = np.pi * b * np.sin(alpha) / laser
  y = 2.*np.pi * ds_abstand * np.sin(alpha) / laser
  return (np.sin(x_skr)/x_skr)**2 * (np.sin(y)/y)**2 - max1[1]
#ds_breite = fsolve(ds_breite_nst, 0.1*es_breite)[0]
ds_breite = 0.1*es_breite # Geht nicht da Spaltabstand zu klein für Messwerte

def I_ds_theo(x):
  alpha = np.arctan2(x,schirmabstand)
  x_skr = np.pi * ds_breite * np.sin(alpha) / laser
  y = 2.*np.pi * ds_abstand * np.sin(alpha) / laser
  return (np.sin(x_skr)/x_skr)**2 * (np.sin(y)/y)**2

x_ds_theo = np.linspace(-4e-4,4e-4,2000)
plt.plot(x_ds_theo,I_ds_theo(x_ds_theo),label='berechnet')
plt.legend()

print()
print(val('DS Breite',ds_breite))

# Berechnung der Beugungsbilder
ds_v = 2
ds_x = np.linspace(-2,2,2000)
es_x = np.linspace(-2,2,2000)

def I_es_theo_skr(x):
  return np.sinc(x)**2

def I_ds_theo_skr(x):
  return np.sinc(x)**2 * np.cos(np.pi*ds_v*x)**2

pltext.initplot(num=5, title='Beugung von DS unter ES', xlabel='pi*x',ylabel='Intensität')
plt.plot(es_x,I_es_theo_skr(es_x),label='ES')
plt.plot(ds_x,I_ds_theo_skr(ds_x),label='DS')
plt.legend()
plt.ylim((0,1.1))

# Berechnung modifizierter Spaltbilder
d=1
g=2*d

def I_es_mod_skr(k):
  return d/np.pi * np.sin(k*d/2)/(k*d/2) * np.cos(es_y*k)

def I_ds_mod_skr(k):
  return d/np.pi*np.cos(k*g/2)*np.sin(k*d/2)/(k*d/2)*np.cos(ds_y*k)

es_Y=np.linspace(-1,1,2000)*d
ds_Y=np.linspace(-2,2,2000)*d

from scipy.integrate import quad

for n in range(1,6):

  es_f_modif = []
  ds_f_modif = []
  for i in range(len(es_Y)):
    es_y = es_Y[i]
    ds_y = ds_Y[i]
    es_res,es_err = quad(I_es_mod_skr,0,2*np.pi*n/d)
    ds_res,ds_err = quad(I_ds_mod_skr,0,2*np.pi*n/d)
    es_f_modif.append(es_res**2)
    ds_f_modif.append(ds_res**2)
  es_f_modif=es_f_modif/np.max(es_f_modif)
  ds_f_modif=ds_f_modif/np.max(ds_f_modif)

  pltext.initplot(num=6+2*(n-1),title='modifizierter ES mit n='+str(n),xlabel='y / d', ylabel='Intensität')
  plt.plot(es_Y,es_f_modif)
  pltext.initplot(num=7+2*(n-1),title='modifizierter DS mit n='+str(n),xlabel='y / d', ylabel='Intensität')
  plt.plot(ds_Y,ds_f_modif)

plt.show()
