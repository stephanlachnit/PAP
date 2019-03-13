# measure version 1.9.1s
from measure import plt,pltext,linreg,sqrt,val,sin,cos,h as h_lit,c,e,deg_to_rad,rad_to_deg,dev,arcsin
from numpy import loadtxt

# Aufgabe 1a
alpha,rate = loadtxt('data/255_1a.txt',unpack=True)
rate_err = sqrt(rate)

pltext.initplot(num=1,title='Abbildung   : Zählrate als Funktion des Winkels',xlabel='Winkel in deg',ylabel='Zählrate')
pltext.plotdata(alpha,rate,rate_err,label='Messwerte',connect=True)
pltext.set_layout(xlim=(2,22.5),ylim=(-0.1e3,1.6e3))

pltext.initplot(num=2,title='Abbildung   : Zählrate als Funktion des Winkels',xlabel='Winkel in deg',ylabel='Zählrate')
pltext.plotdata(alpha,rate,rate_err,label='Messwerte')
[slope,slope_err,yitc,yitc_err] = linreg(alpha[12:18],rate[12:18],rate_err[12:18],plot=True,prange=(4,7))
pltext.set_layout(xlim=(4,7),ylim=(-100,400))

a1 = -yitc / slope * deg_to_rad
a1_err = 1 / slope * sqrt(yitc_err**2 + (yitc / slope * slope_err)**2) * deg_to_rad

d = 201.4e-12
U = 35e3

l = 2 * d * sin(a1)
l_err = 2 * d * cos(a1) * a1_err

h = l * e * U / c
h_err = l_err * e * U / c

a2 = arcsin(l/d)
a2_err = 1/(d * sqrt(1 - (l/d)**2)) * l_err

print('\nAufgabe 1a:\n')
print(val(a1*rad_to_deg,a1_err*rad_to_deg,'a1'))
print(val(l,l_err,'l'))
print(val(h,h_err,'h'))
print(dev(h,h_err,h_lit,name='Abw',perc=True))
print(val(a2*rad_to_deg,a2_err*rad_to_deg,'a2'))

# Aufgabe 1b
alpha,rate = loadtxt('data/255_1b_1o.txt',unpack=True)
rate_err = sqrt(rate)

pltext.initplot(num=3)
pltext.plotdata(alpha,rate,rate_err)
pltext.set_layout()

# Plot
print()
plt.show()
