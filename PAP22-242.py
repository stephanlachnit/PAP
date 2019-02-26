# measure version 1.8.8s
from measure import npfarray,plt,pltext,sqrt,mv,dtot_mv,tbl,lst,sig,dev

# Aufgabe 1
U1g = npfarray([-0.2,-0.15,-0.10,-0.05,0.04,0.08,0.12,0.18])
U1g_dsys = npfarray([3,1,1,1,1,1,1,1])*1e-3
Uag_48k7 = npfarray([3.32,2.60,1.80,1.00,-0.4,-1.0,-1.7,-2.7])
Uag_48k7_dsys = npfarray([50,10,10,10,1,10,100,100])*1e-3
Uag_274k = npfarray([14.6,14.4,9.6,5.0,-3.2,-6.90,-10.8,-12.8])
Uag_274k_dsys = npfarray([0.1,0.1,0.01,0.01,0.01,0.1,0.1,0.1])

U1w = npfarray([1.0,0.852,0.748,0.500,0.300,0.150]) / 10
U1w_dsys = npfarray([10,5,5,5,5,3])*1e-3 / 10
Uaw_274k = npfarray([8.76,7.48,6.46,4.38,2.61,1.28])
Uaw_274k_dsys = npfarray([0.05,0.05,0.03,0.03,0.02,0.02])
Uaw_680k = npfarray([21.9,18.7,16.1,10.8,6.52,3.25])
Uaw_680k_dsys = npfarray([0.1,0.1,0.1,0.1,0.05,0.02])

V_eg_48k7 = -1 * Uag_48k7 / U1g
V_eg_48k7_dsys = 1/abs(U1g) * sqrt(Uag_48k7_dsys**2 + (Uag_48k7 * U1g_dsys / U1g)**2)
V_eg_48k7_mv = mv(V_eg_48k7)
V_eg_48k7_mv_dtot = dtot_mv(V_eg_48k7,V_eg_48k7_dsys)

V_eg_274k = -1 * Uag_274k / U1g
V_eg_274k_dsys = 1/abs(U1g) * sqrt(Uag_274k_dsys**2 + (Uag_274k * U1g_dsys / U1g)**2)
V_eg_274k_mv = mv(V_eg_274k)
V_eg_274k_mv_dtot = dtot_mv(V_eg_274k,V_eg_274k_dsys)

V_ew_274k = Uaw_274k / U1w
V_ew_274k_dsys = 1/U1w * sqrt(Uaw_274k_dsys**2 + (Uaw_274k * U1w_dsys / U1w)**2)
V_ew_274k_mv = mv(V_ew_274k)
V_ew_274k_mv_dtot = dtot_mv(V_ew_274k,V_ew_274k_dsys)

V_ew_680k = Uaw_680k / U1w
V_ew_680k_dsys = 1/U1w * sqrt(Uaw_680k_dsys**2 + (Uaw_680k * U1w_dsys / U1w)**2)
V_ew_680k_mv = mv(V_ew_680k)
V_ew_680k_mv_dtot = dtot_mv(V_ew_680k,V_ew_680k_dsys)

R_E = 3e3
R_G = npfarray([48.7e3,274e3,680e3])
V_t = R_G / R_E

print('\nAufgabe 1:\n')
print(tbl([['Widerstand','48k7 (g)','274k (g)','274k (w)','680k (w)'],lst([V_t[0],V_t[1],V_t[1],V_t[2]],name='V_t'),lst([V_eg_48k7_mv,V_eg_274k_mv,V_ew_274k_mv,V_ew_680k_mv],[V_eg_48k7_mv_dtot,V_eg_274k_mv_dtot,V_ew_274k_mv_dtot,V_ew_680k_mv_dtot],'V_e'),['Abw',sig('',V_eg_48k7_mv,V_eg_48k7_mv_dtot,V_t[0],perc=True),sig('',V_eg_274k_mv,V_eg_274k_mv_dtot,V_t[1],perc=True),sig('',V_ew_274k_mv,V_ew_274k_mv_dtot,V_t[1],perc=True),sig('',V_ew_680k_mv,V_ew_680k_mv_dtot,V_t[2],perc=True)]]))
print('\nEinzelwerte für 48k7 (g):\n')
print(tbl([lst(V_eg_48k7,V_eg_48k7_dsys,'Verstärkung'),dev(V_eg_48k7,V_eg_48k7_dsys,V_t[0],name='Abw',perc=True)]))
print('\nEinzelwerte für 274k (g):\n')
print(tbl([lst(V_eg_274k,V_eg_274k_dsys,'Verstärkung'),dev(V_eg_274k,V_eg_274k_dsys,V_t[1],name='Abw',perc=True)]))
print('\nEinzelwerte für 274k (w):\n')
print(tbl([lst(V_ew_274k,V_ew_274k_dsys,'Verstärkung'),dev(V_ew_274k,V_ew_274k_dsys,V_t[1],name='Abw',perc=True)]))
print('\nEinzelwerte für 680k (w):\n')
print(tbl([lst(V_ew_680k,V_ew_680k_dsys,'Verstärkung'),dev(V_ew_680k,V_ew_680k_dsys,V_t[2],name='Abw',perc=True)]))

# Plots
f_uncert = 50e-6

f_1 = npfarray([0.3,0.6,0.9,3,6,9,30,60,90,150,200,300])*1e3
f_1_dsys = f_1 * f_uncert
U_A1 = npfarray([6.76,6.70,6.68,5.46,3.82,2.78,0.900,0.456,0.306,0.187,0.140,0.095])
U_A1_dsys = npfarray([0.01,0.03,0.05,0.03,0.03,0.05,0.005,0.005,0.003,0.003,0.001,0.001])
U_A2 = npfarray([2.65,2.64,2.67,2.56,2.30,2.00,0.852,0.450,0.304,0.186,0.139,0.094])
U_A2_dsys = npfarray([0.02,0.01,0.02,0.05,0.01,0.03,0.003,0.003,0.003,0.003,0.003,0.001])
U_A3 = npfarray([1.52,1.52,1.52,1.54,1.53,1.52,1.36,1.06,0.836,0.554,0.428,0.293])
U_A3_dsys = npfarray([0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.005,0.003,0.005,0.003])
U_A4 = npfarray([1.53,1.53,1.51,1.39,1.11,0.884,0.316,0.164,0.109,0.0664,0.0508,0.0348])
U_A4_dsys = npfarray([0.01,0.01,0.04,0.02,0.01,0.004,0.002,0.002,0.001,0.0004,0.0004,0.0002])

f_2 = npfarray([0.3,0.6,0.9,1.5,3.0,6.0,9.0,11.5,14.0,16.0,18.0,20.0])*1e3
f_2_dsys = f_2 * f_uncert
U_A5 = npfarray([424,744,984,1250,1460,1540,1530,1520,1500,1500,1480,1460])*1e-3
U_A5_dsys = npfarray([1,1,1,10,10,10,10,10,10,10,10,10])*1e-3

Vss = [0.3,1.0]
Vss_dsys = [5e-3,20e-3]
V_A1 = U_A1 / Vss[0]
V_A1_dsys = 1/Vss[0] * sqrt(U_A1_dsys**2 + (U_A1 * Vss_dsys[0] / Vss[0])**2)
V_A2 = U_A2 / Vss[0]
V_A2_dsys = 1/Vss[0] * sqrt(U_A2_dsys**2 + (U_A2 * Vss_dsys[0] / Vss[0])**2)
V_A3 = U_A3 / Vss[1]
V_A3_dsys = 1/Vss[1] * sqrt(U_A3_dsys**2 + (U_A3 * Vss_dsys[1] / Vss[1])**2)
V_A4 = U_A4 / Vss[1]
V_A4_dsys = 1/Vss[1] * sqrt(U_A4_dsys**2 + (U_A4 * Vss_dsys[1] / Vss[1])**2)
V_A5 = U_A5 / Vss[1]
V_A5_dsys = 1/Vss[1] * sqrt(U_A5_dsys**2 + (U_A5 * Vss_dsys[1] / Vss[1])**2)

pltext.initplot(num=1,title='Abbildung   : Frequenzgang des Verstärkers (Spannung)',xlabel='Frequenz in Hz',ylabel='Spannung in V',scale='loglog')
pltext.plotdata(f_1,U_A1,U_A1_dsys,f_1_dsys,label='680k',caps=False)
pltext.plotdata(f_1,U_A2,U_A2_dsys,f_1_dsys,label='274k',caps=False)
pltext.plotdata(f_1,U_A3,U_A3_dsys,f_1_dsys,label='48k7',caps=False)
pltext.plotdata(f_1,U_A4,U_A4_dsys,f_1_dsys,label='48k7 mit 560pF Parallelkapazität',caps=False)
pltext.plotdata(f_2,U_A5,U_A5_dsys,f_2_dsys,label='48k7 mit 47nF Eingangskapazität',caps=False)
pltext.set_layout(legend=True,xlim=(260,3.4e5),ylim=(3e-2,8))

pltext.initplot(num=2,title='Abbildung   : Frequenzgang des Verstärkers (Verstärkung)',xlabel='Frequenz in Hz',ylabel='Verstärkung',scale='loglog')
pltext.plotdata(f_1,V_A1,V_A1_dsys,f_1_dsys,label='680k',caps=False)
pltext.plotdata(f_1,V_A2,V_A2_dsys,f_1_dsys,label='274k',caps=False)
pltext.plotdata(f_1,V_A3,V_A3_dsys,f_1_dsys,label='48k7',caps=False)
pltext.plotdata(f_1,V_A4,V_A4_dsys,f_1_dsys,label='48k7 mit 560pF Parallelkapazität',caps=False)
pltext.plotdata(f_2,V_A5,V_A5_dsys,f_2_dsys,label='48k7 mit 47nF Eingangskapazität',caps=False)
pltext.set_layout(legend=True,xlim=(260,3.4e5),ylim=(3e-2,3e2))

print()
plt.show()
