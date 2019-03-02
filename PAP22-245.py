# measure version 1.8.9s
from measure import sqrt,val,npfarray,pltext,plt,linreg,pi,dev

r_h = 0.295 / 2
s_h = 0.147
N_h = 124

N_i = 4000
A_i = 41.7e-4

# Aufgabe 2
I_f = 4.00
I_f_dsys = sqrt(0.02**2 + 0.005**2 + (0.012 * I_f)**2)
f_f = npfarray([3.1,6.1,8.8,11.8,14.9])
f_f_dsys = npfarray([0.02,0.10,0.10,0.10,0.10])
Vss_f = npfarray([0.88,2.70,4.80,7.00,9.30])
Vss_f_dsys = npfarray([0.02,0.05,0.10,0.10,0.10])
Vind_f = Vss_f / 2
Vind_f_dsys = Vss_f_dsys / 2

f_I = 9.46
f_I_dsys = 0.10
I_I = npfarray([0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5])
I_I_dsys = sqrt(npfarray([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])**2 + 0.005**2 + (0.012 * I_I)**2)
Vss_I = npfarray([0.724,1.41,2.04,2.68,3.32,4.00,4.65,5.32,6.12])
Vss_I_dsys = npfarray([0.01,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
Vind_I = Vss_I / 2
Vind_I_dsys = Vss_I_dsys / 2

pltext.initplot(num=1,title='Abbildung   : Induktionsspannung als Funktion der Frequenz',xlabel='Drehfrequenz in Hz',ylabel='Induktionsspannung in V')
[slope,dslope] = linreg(f_f,Vind_f,Vind_f_dsys,f_f_dsys,plot=True)[0:2]
pltext.set_layout(legend=True,xlim=(2,16),ylim=(0,5))

B_f = slope / (2*pi * A_i * N_i)
B_f_dsys = dslope / (2*pi * A_i * N_i)

mu_0 = 4*pi * 1e-7
B_theo = (8 / sqrt(125)) * mu_0 * I_f * N_h / s_h
B_theo_dsys = (8 / sqrt(125)) * mu_0 * I_f_dsys * N_h / s_h

pltext.initplot(num=2,title='Abbildung   : Induktionsspannung als Funktion des Stroms',xlabel='Spulenstrom in A',ylabel='Induktionsspannung in V')
linreg(I_I,Vind_I,Vind_I_dsys,I_I_dsys,plot=True)
pltext.set_layout(legend=True,xlim=(0,5),ylim=(0,3.5))

print('\nAufgabe 2:\n')
print(val(B_f,B_f_dsys,'B_exp '))
print(val(B_theo,B_theo_dsys,'B_theo'))
print(dev(B_f,B_f_dsys,B_theo,B_theo_dsys,'Abw'))

# Aufgabe 3
w = npfarray([0,30,60,90])
w_dsys = npfarray([2.5,2.5,2.5,2.5])
Vss_w = npfarray([1.48,1.29,0.73,0.03])
Vss_w_dsys = npfarray([0.03,0.02,0.01,0.02])
Vind_w =  Vss_w / 2
Vind_w_dsys = Vss_w_dsys / 2

pltext.initplot(num=3,title='Abbildung   : Induktionsspannung als Funktion des Winkel',xlabel='Winkel in deg',ylabel='Induktionsspannung in V')
pltext.plotdata(w,Vind_w,Vind_w_dsys,w_dsys)
pltext.set_layout()

# Plots
print()

plt.show()
