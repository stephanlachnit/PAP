# measure version 1.9.2s
from measure import plt,pltext,linreg,sqrt,val,sin,cos,h as h_lit,NA as NA_lit,c,e,deg_to_rad,rad_to_deg,dev,arcsin,spcurvefit,exp,pi,nplinspace,ln,tbl,lst,npfarray,mv,dtot_mv,loadtxt

# Aufgabe 1a
alpha,rate = loadtxt('data/255_1a.txt',unpack=True)
rate_err = sqrt(rate)

pltext.initplot(num=1,title='Abbildung   : Zählrate als Funktion des Winkels (LiF)',xlabel='Winkel in deg',ylabel='Zählrate in 1/s')
pltext.plotdata(alpha,rate,rate_err,label='Messwerte',connect=True)
pltext.set_layout(xlim=(2,22.5),ylim=(-0.1e3,1.6e3))

pltext.initplot(num=2,title='Abbildung   : Zählrate als Funktion des Winkels (LiF)',xlabel='Winkel in deg',ylabel='Zählrate in 1/s')
pltext.plotdata(alpha,rate,rate_err,label='Messwerte')
[slope,slope_err,yitc,yitc_err] = linreg(alpha[12:18],rate[12:18],rate_err[12:18],plot=True,prange=(4,7))
pltext.set_layout(xlim=(4,7),ylim=(-100,400))

a1_uncorr = -yitc / slope * deg_to_rad
a1_uncorr_err = 1 / slope * sqrt(yitc_err**2 + (yitc / slope * slope_err)**2) * deg_to_rad

Ug = mv(rate[0:12])
Ug_err = dtot_mv(rate[0:12],rate_err[0:12])

a1 = (1 - Ug / yitc) * a1_uncorr
a1_err = sqrt((a1_uncorr / yitc)**2 * (Ug_err**2 + (Ug / yitc * yitc)**2) + ((1- Ug / yitc) * a1_uncorr_err)**2)

d = 201.4e-12
U = 35e3

l = 2 * d * sin(a1)
l_err = 2 * d * cos(a1) * a1_err

h = l * e * U / c
h_err = l_err * e * U / c

a2 = arcsin(l/d)
a2_err = 1/(d * sqrt(1 - (l/d)**2)) * l_err

print('\nAufgabe 1a:\n')
print(val(a1_uncorr*rad_to_deg,a1_uncorr_err*rad_to_deg,'a1 (uncorr)'))
print(val(Ug,Ug_err,'Ug'))
print(val(a1*rad_to_deg,a1_err*rad_to_deg,'a1'))
print(val(l,l_err,'l'))
print(val(h,h_err,'h'))
print(dev(h,h_err,h_lit,name='Abw',perc=True))
print(val(a2*rad_to_deg,a2_err*rad_to_deg,'a2'))

# Aufgabe 1b
alpha_1o,rate_1o = loadtxt('data/255_1b_1o.txt',unpack=True)
rate_1o_err = sqrt(rate_1o)
alpha_2o,rate_2o = loadtxt('data/255_1b_2o.txt',unpack=True)
rate_2o_err = sqrt(rate_2o)

def gauss(x,A,mu,sig,Ug):
  return A / (sqrt(2*pi) * sig) * exp(-(x-mu)**2 / (2 * sig**2)) + Ug

l_kb_lit = 63.1e-12
l_ka_lit = 71.1e-12

p_opt_kb1o,p_err_kb1o = spcurvefit(gauss,alpha_1o[5:14],rate_1o[5:14],rate_1o_err[5:14],[100,9,0.1,300])
p_opt_ka1o,p_err_ka1o = spcurvefit(gauss,alpha_1o[-15:-5],rate_1o[-15:-5],rate_1o_err[-15:-5],[400,10,0.1,200])

fwhm_kb1o = 2*sqrt(2*ln(2)) * p_opt_kb1o[2]
fwhm_kb1o_err = 2*sqrt(2*ln(2)) * p_err_kb1o[2]
fwhm_ka1o = 2*sqrt(2*ln(2)) * p_opt_ka1o[2]
fwhm_ka1o_err = 2*sqrt(2*ln(2)) * p_err_ka1o[2]

l_kb1o = 2 * d * sin(p_opt_kb1o[1]*deg_to_rad)
l_kb1o_err = 2 * d * cos(p_opt_kb1o[1]*deg_to_rad) * p_err_kb1o[1] * deg_to_rad
l_ka1o = 2 * d * sin(p_opt_ka1o[1]*deg_to_rad)
l_ka1o_err = 2 * d * cos(p_opt_ka1o[1]*deg_to_rad) * p_err_ka1o[1] * deg_to_rad

p_opt_kb2o,p_err_kb2o = spcurvefit(gauss,alpha_2o[2:15],rate_2o[2:15],rate_2o_err[2:15],[20,18.3,0.1,60])
p_opt_ka2o,p_err_ka2o = spcurvefit(gauss,alpha_2o[-15:-3],rate_2o[-15:-3],rate_2o_err[-15:-3],[60,20.7,0.1,50])

fwhm_kb2o = 2*sqrt(2*ln(2)) * p_opt_kb2o[2]
fwhm_kb2o_err = 2*sqrt(2*ln(2)) * p_err_kb2o[2]
fwhm_ka2o = 2*sqrt(2*ln(2)) * p_opt_ka2o[2]
fwhm_ka2o_err = 2*sqrt(2*ln(2)) * p_err_ka2o[2]

l_kb2o = d * sin(p_opt_kb2o[1]*deg_to_rad)
l_kb2o_err = d * cos(p_opt_kb2o[1]*deg_to_rad) * p_err_kb2o[1] * deg_to_rad
l_ka2o = d * sin(p_opt_ka2o[1]*deg_to_rad)
l_ka2o_err = d * cos(p_opt_ka2o[1]*deg_to_rad) * p_err_ka2o[1] * deg_to_rad

x_kb1o_array = nplinspace(8.5,9.4)
x_ka1o_array = nplinspace(9.6,10.7)
pltext.initplot(num=3,title='Abbildung   : Extrema erster Ordnung (LiF)',xlabel='Winkel in deg',ylabel='Zählrate in 1/s')
pltext.plotdata(alpha_1o,rate_1o,rate_1o_err,label='Messwerte')
plt.plot(x_kb1o_array,gauss(x_kb1o_array,*p_opt_kb1o),label=r'$K_\beta$ Fit')
plt.plot(x_ka1o_array,gauss(x_ka1o_array,*p_opt_ka1o),label=r'$K_\alpha$ Fit')
pltext.set_layout(xlim=(7.9,11.1),ylim=(150,1600))

x_kb2o_array = nplinspace(17.7,18.8)
x_ka2o_array = nplinspace(20.1,21.2)
pltext.initplot(num=4,title='Abbildung   : Extrema zweiter Ordnung (LiF)',xlabel='Winkel in deg',ylabel='Zählrate in 1/s')
pltext.plotdata(alpha_2o,rate_2o,rate_2o_err,label='Messwerte')
plt.plot(x_kb2o_array,gauss(x_kb2o_array,*p_opt_kb2o),label=r'$K_\beta$ Fit')
plt.plot(x_ka2o_array,gauss(x_ka2o_array,*p_opt_ka2o),label=r'$K_\alpha$ Fit')
pltext.set_layout(xlim=(17.4,21.6))

print('\nAufgabe 1b\n')
print(tbl([['Peak:',' A',' mu',' sig',' Ug',' FWHM',' l',' Abw'],lst([*p_opt_kb1o,fwhm_kb1o,l_kb1o],[*p_err_kb1o,fwhm_kb1o_err,l_kb1o_err],'kb1o')+[dev(l_kb1o,l_kb1o_err,l_kb_lit,perc=True)],lst([*p_opt_ka1o,fwhm_ka1o,l_ka1o],[*p_err_ka1o,fwhm_ka1o_err,l_ka1o_err],'ka1o')+[dev(l_ka1o,l_ka1o_err,l_ka_lit,perc=True)],lst([*p_opt_kb2o,fwhm_kb2o,l_kb2o],[*p_err_kb2o,fwhm_kb2o_err,l_kb2o_err],'kb2o')+[dev(l_kb2o,l_kb2o_err,l_kb_lit,perc=True)],lst([*p_opt_ka2o,fwhm_ka2o,l_ka2o],[*p_err_ka2o,fwhm_ka2o_err,l_ka2o_err],'ka2o')+[dev(l_ka2o,l_ka2o_err,l_ka_lit,perc=True)]]))

# Aufgabe 1c
U = npfarray([20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35])
rate = npfarray([1.00,1.90,2.25,8.00,37.55,72.20,104.0,134.4,160.2,180.5,212.5,245.6,264.6,284.8,314.1,339.1])
rate_err = sqrt(rate)

pltext.initplot(num=5,title='Abbildung   : Zählrate als Funktion der Spannung (LiF)',xlabel='Spannung in kV',ylabel='Zählrate in 1/s')
pltext.plotdata(U,rate,rate_err,label='Messwerte')
[slope,slope_err,yitc,yitc_err] = linreg(U[3:],rate[3:],rate_err[3:],plot=True,prange=(19.5,35.5))
pltext.set_layout(xlim=(19.5,35.5),ylim=(-25,375))

U0 = -yitc / slope * 1e3
U0_err = 1 / slope * sqrt(yitc_err**2 + (yitc / slope * slope_err)**2) * 1e3

l = 2 * d * sin(7.5 * deg_to_rad)
h = l * e * U0 / c
h_err = l * e / c * U0_err

print('\nAufgabe 1c\n')
print(val(U0,U0_err,'U0'))
print(val(h,h_err,'h'))
print(dev(h,h_err,h_lit,name='Abw',perc=True))

# Aufgabe 2
alpha,rate = loadtxt('data/255_2.txt',unpack=True)
rate_err = sqrt(rate)

p_opt_kb1o,p_err_kb1o = spcurvefit(gauss,alpha[14:19],rate[14:19],rate_err[14:19],[100,6.2,0.1,400])
p_opt_ka1o,p_err_ka1o = spcurvefit(gauss,alpha[18:24],rate[18:24],rate_err[18:24],[300,7.0,0.1,300])
p_opt_kb2o,p_err_kb2o = spcurvefit(gauss,alpha[-30:-23],rate[-30:-23],rate_err[-30:-23],[40,12.7,0.1,80])
p_opt_ka2o,p_err_ka2o = spcurvefit(gauss,alpha[-23:-12],rate[-23:-12],rate_err[-23:-12],[100,14.4,0.1,70])

g_kb1o = l_kb_lit / (2 * sin(p_opt_kb1o[1]*deg_to_rad))
g_kb1o_err = l_kb_lit / (2 * sin(p_opt_kb1o[1]*deg_to_rad)**2) * cos(p_opt_kb1o[1]*deg_to_rad) * p_err_kb1o[1]*deg_to_rad
g_ka1o = l_ka_lit / (2 * sin(p_opt_ka1o[1]*deg_to_rad))
g_ka1o_err = l_ka_lit / (2 * sin(p_opt_ka1o[1]*deg_to_rad)**2) * cos(p_opt_ka1o[1]*deg_to_rad) * p_err_ka1o[1]*deg_to_rad
g_kb2o = l_kb_lit / sin(p_opt_kb2o[1]*deg_to_rad)
g_kb2o_err = l_kb_lit / sin(p_opt_kb2o[1]*deg_to_rad)**2 * cos(p_opt_kb2o[1]*deg_to_rad) * p_err_kb2o[1]*deg_to_rad
g_ka2o = l_ka_lit / sin(p_opt_ka2o[1]*deg_to_rad)
g_ka2o_err = l_ka_lit / sin(p_opt_ka2o[1]*deg_to_rad)**2 * cos(p_opt_ka2o[1]*deg_to_rad) * p_err_ka2o[1]*deg_to_rad

g = mv([g_kb1o,g_ka1o,g_kb2o,g_ka2o])
g_err = dtot_mv([g_kb1o,g_ka1o,g_kb2o,g_ka2o],[g_kb1o_err,g_ka1o_err,g_kb2o_err,g_ka2o_err])

M = 58.44e-3
roh = 2.164e3

NA = 1/2 * M / (roh * g**3)
NA_err = 3/2 * M / (roh * g**4) * g_err

x_kb1o_array = nplinspace(5.8,6.6)
x_ka1o_array = nplinspace(6.6,7.6)
x_kb2o_array = nplinspace(12.2,13.4)
x_ka2o_array = nplinspace(13.6,15.6)
pltext.initplot(num=6,title='Abbildung   : Zählrate als Funktion des Winkel (NaCl)',xlabel='Winkel in deg',ylabel='Zählrate in 1/s')
pltext.plotdata(alpha,rate,rate_err,label='Messwerte')
plt.plot(x_kb1o_array,gauss(x_kb1o_array,*p_opt_kb1o),label=r'$K_\beta$ Fit (1. Ordnung)')
plt.plot(x_ka1o_array,gauss(x_ka1o_array,*p_opt_ka1o),label=r'$K_\alpha$ Fit (1. Ordnung)')
plt.plot(x_kb2o_array,gauss(x_kb2o_array,*p_opt_kb2o),label=r'$K_\beta$ Fit (2. Ordnung)')
plt.plot(x_ka2o_array,gauss(x_ka2o_array,*p_opt_ka2o),label=r'$K_\alpha$ Fit (2. Ordnung)')
pltext.set_layout(xlim=(2.75,18.25),ylim=(-40,1640))

print('\nAufgabe 2:\n')
print(tbl([['Peak:',' A',' mu',' sig',' Ug',' g'],lst([*p_opt_kb1o,g_kb1o],[*p_err_kb1o,g_kb1o_err],'kb1o'),lst([*p_opt_ka1o,g_ka1o],[*p_err_ka1o,g_ka1o_err],'ka1o'),lst([*p_opt_kb2o,g_kb2o],[*p_err_kb2o,g_kb2o_err],'kb2o'),lst([*p_opt_ka2o,g_ka2o],[*p_err_ka2o,g_ka2o_err],'ka2o')]))
print()
print(val(g,g_err,'g'))
print(val(NA,NA_err,'NA'))
print(dev(NA,NA_err,NA_lit,name='Abw',perc=True))

# Plot
print()
plt.show()
