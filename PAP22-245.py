# measure version 1.8.10s
from measure import sqrt,val,npfarray,pltext,plt,linreg,pi,dev,arccos

r_h = 0.295 / 2
s_h = 0.147
N_h = 124

N_i = 4000
A_i = 41.7e-4

# Aufgabe 2a
I_2a = 4.00
I_2a_dsys = sqrt(0.02**2 + 0.005**2 + (0.012 * I_2a)**2)
f_2a = npfarray([3.1,6.1,8.8,11.8,14.9])
f_2a_dsys = npfarray([0.02,0.10,0.10,0.10,0.10])
Vss_2a = npfarray([0.88,2.70,4.80,7.00,9.30])
Vss_2a_dsys = npfarray([0.02,0.05,0.10,0.10,0.10])

Uind_2a = Vss_2a / 2
Uind_2a_dsys = Vss_2a_dsys / 2

pltext.initplot(num=1,title='Abbildung   : Induktionsspannung als Funktion der Frequenz',xlabel='Drehfrequenz in Hz',ylabel='Induktionsspannung in V')
[slope_2a,dslope_2a] = linreg(f_2a,Uind_2a,Uind_2a_dsys,f_2a_dsys,plot=True,prange=(2,16))[0:2]
pltext.set_layout(legend=True,xlim=(2,16),ylim=(0,5))

B_exp_2a = slope_2a / (2*pi * A_i * N_i)
B_exp_2a_dsys = dslope_2a / (2*pi * A_i * N_i)

mu_0 = 4*pi * 1e-7
B_theo_2a = (8 / sqrt(125)) * mu_0 * I_2a * N_h / s_h
B_theo_2a_dsys = (8 / sqrt(125)) * mu_0 * I_2a_dsys * N_h / s_h

print('\nAufgabe 2a:\n')
print(val(B_exp_2a,B_exp_2a_dsys,'B_exp '))
print(val(B_theo_2a,B_theo_2a_dsys,'B_theo'))
print(dev(B_exp_2a,B_exp_2a_dsys,B_theo_2a,B_theo_2a_dsys,'Abw'))

# Aufgabe 2b
f_2b = 9.46
f_2b_dsys = 0.10
I_2b = npfarray([0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5])
I_2b_dsys = sqrt(npfarray([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])**2 + 0.005**2 + (0.012 * I_2b)**2)
Vss_2b = npfarray([0.724,1.41,2.04,2.68,3.32,4.00,4.65,5.32,6.12])
Vss_2b_dsys = npfarray([0.01,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])

Uind_2b = Vss_2b / 2
Uind_2b_dsys = Vss_2b_dsys / 2

pltext.initplot(num=2,title='Abbildung   : Induktionsspannung als Funktion des Stroms',xlabel='Spulenstrom in A',ylabel='Induktionsspannung in V')
linreg(I_2b,Uind_2b,Uind_2b_dsys,I_2b_dsys,plot=True,prange=(0,5))
pltext.set_layout(legend=True,xlim=(0,5),ylim=(0,3.5))

# Aufgabe 3a
w_3a = npfarray([0,30,60,90])
w_3a_dsys = npfarray([2.5,2.5,2.5,2.5])
Vss_3a = npfarray([1.48,1.29,0.73,0.03])
Vss_3a_dsys = npfarray([0.03,0.02,0.01,0.02])

Uind_3a =  Vss_3a / 2
Uind_3a_dsys = Vss_3a_dsys / 2

pltext.initplot(num=3,title='Abbildung   : Induktionsspannung als Funktion des Winkel',xlabel='Winkel in deg',ylabel='Induktionsspannung in V')
pltext.plotdata(w_3a,Uind_3a,Uind_3a_dsys,w_3a_dsys,connect=True)
pltext.set_layout(xlim=(-10,100),ylim=(0,0.8))

# Aufgabe 3b
f_3b = npfarray([20.3,40.4,60.2,80.1,100.2,120.0,142.3,165.0,180.3,200.5,404.5,595.5,802.5,1006,1206,1404,1603,1784,2025])
f_3b_dsys = npfarray([0.1,0.1,0.2,0.1,0.3,0.2,0.2,0.3,0.3,0.2,0.3,0.4,1.5,2,2,2,3,4,3])
Vss_ind_3b = npfarray([1.14,2.05,2.65,3.02,3.26,3.42,3.52,3.59,3.64,3.66,3.82,3.84,3.86,3.86,3.92,3.96,4.00,4.06,4.12])
Vss_ind_3b_dsys = npfarray([0.05,0.05,0.05,0.05,0.02,0.05,0.02,0.05,0.05,0.02,0.05,0.03,0.02,0.02,0.05,0.05,0.05,0.05,0.02])
I_3b = npfarray([244.0,217.4,187.1,160.1,137.9,120.4,104.9,92.6,83.3,77.7,39.9,27.41,20.42,16.32,13.64,11.72,10.28,9.23,8.13])*1e-3
I_3b_dsys = sqrt((npfarray([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])*1e-3)**2 + (0.03e-3)**2 + (0.01 * I_3b)**2)
Vss_h_3b = npfarray([5.50,5.65,5.71,5.72,5.70,5.68,5.68,5.68,5.70,5.70,5.72,5.72,5.72,5.72,5.72,5.74,5.72,5.78,5.76])
Vss_h_3b_dsys = npfarray([0.10,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])

Uind_3b = Vss_ind_3b / 2
Uind_3b_dsys = Vss_ind_3b_dsys / 2
Uh_3b = Vss_h_3b / 2
Uh_3b_dsys = Vss_h_3b_dsys / 2

Uind_Uh_ratio = Uind_3b / Uh_3b
Uind_Uh_ratio_dsys = 1 / Uh_3b * sqrt(Uind_3b_dsys**2 + (Uind_3b * Uh_3b_dsys / Uh_3b)**2)

pltext.initplot(num=4,title=r'Abbildung   : Verhältnis $U_{ind}$ / $U_h$ als Funktion der Frequenz',xlabel='Frequenz in Hz',ylabel=r'$U_{ind}$ / $U_h$')
pltext.plotdata(f_3b,Uind_Uh_ratio,Uind_Uh_ratio_dsys,f_3b_dsys)
pltext.set_layout(xlim=(0,2.25e3),ylim=(0.1,0.8))

# Aufgabe 3c
R_3c = Uh_3b / I_3b
R_3c_dsys = 1 / I_3b * sqrt(Uh_3b_dsys**2 + (Uh_3b * I_3b_dsys / I_3b)**2)

pltext.initplot(num=5,title='Abbildung   : Spulenwiderstand als Funktion der Frequenz',xlabel='Frequenz in Hz',ylabel='Spulenwiderstand in Ohm')
[slope_3c,dslope_3c] = linreg(f_3b,R_3c,R_3c_dsys,f_3b_dsys,plot=True,prange=(0,2.25e3),frange=range(4,19))[0:2]
pltext.set_layout(legend=True,xlim=(0,2.25e3),ylim=(0,400))

L_3c = slope_3c / (2*pi)
L_3c_dsys = dslope_3c / (2*pi)

print('\nAufgabe 3c:\n')
print(val(L_3c,L_3c_dsys,'L'))

# Aufgabe 4a
f_4a = 14.7
f_4a_dsys = 0.1
Uind_4a = 112e-3 / 2
Uind_4a_dsys = 2e-3 / 2

B_4a = Uind_4a / (2*pi * N_i * A_i * f_4a)
B_4a_dsys = 1 / (2*pi * N_i * A_i * f_4a) * sqrt(Uind_4a_dsys**2 + (Uind_4a * f_4a_dsys / f_4a)**2)

print('\nAufgabe 4a:\n')
print(val(B_4a,B_4a_dsys,'B'))

# Aufgab 4b
f_4b = 14.3
f_4b_dsys = 0.1
I_4b = 48e-3
I_4b_dsys = sqrt((0.1e-3)**2 + (0.03e-3)**2 + (0.01 * I_4b)**2)
Uind_4b = 23.2e-3 / 2
Uind_4b_dsys = 0.5e-3 / 2

B_4b = Uind_4b / (2*pi * N_i * A_i * f_4b)
B_4b_dsys = 1 / (2*pi * N_i * A_i * f_4b) * sqrt(Uind_4b_dsys**2 + (Uind_4b * f_4b_dsys / f_4b)**2)

rad_to_deg = 360 / (2*pi)
w_4b = arccos(B_4b / B_4a) * rad_to_deg
w_4b_dsys = 1 / (B_4a * sqrt(1 - (B_4b / B_4a)**2)) * sqrt(B_4b_dsys**2 + (B_4b * B_4a_dsys / B_4a)**2) * rad_to_deg
w_4b_theo = 66

print('\nAufgabe 4b:\n')
print(val(B_4b,B_4b_dsys,'B'))
print(val(w_4b,w_4b_dsys,'w'))
print(dev(w_4b,w_4b_dsys,w_4b_theo,name='Abw',perc=True))

# Plots
print()
plt.show()
