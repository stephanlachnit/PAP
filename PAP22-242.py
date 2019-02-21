# measure version 1.8.7
from measure import npfarray,plt,pltext

f_1 = npfarray([0.3,0.6,0.9,3,6,9,30,60,90,150,200,300])*1e3
U_A1 = npfarray([6.76,6.70,6.68,5.46,3.82,2.78,0.900,0.456,0.306,0.187,0.140,0.095])
U_A1_dsys = npfarray([0.01,0.03,0.05,0.03,0.03,0.05,0.005,0.005,0.003,0.003,0.001,0.001])
U_A2 = npfarray([2.65,2.64,2.67,2.56,2.30,2.00,0.852,0.450,0.304,0.186,0.139,0.094])
U_A2_dsys = npfarray([0.02,0.01,0.02,0.05,0.01,0.03,0.003,0.003,0.003,0.003,0.003,0.001])
U_A3 = npfarray([1.52,1.52,1.52,1.54,1.53,1.52,1.36,1.06,0.836,0.554,0.428,0.293])
U_A3_dsys = npfarray([0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.005,0.003,0.005,0.003])
U_A4 = npfarray([1.53,1.53,1.51,1.39,1.11,0.884,0.316,0.164,0.109,0.0664,0.0508,0.0348])
U_A4_dsys = npfarray([0.01,0.01,0.04,0.02,0.01,0.004,0.002,0.002,0.001,0.0004,0.0004,0.0002])

f_2 = npfarray([0.3,0.6,0.9,1.5,3.0,6.0,9.0,11.5,14.0,16.0,18.0,20.0])*1e3
U_A5 = npfarray([424,744,984,1250,1460,1540,1530,1520,1500,1500,1480,1460])*1e-3
U_A5_dsys = npfarray([1,1,1,10,10,10,10,10,10,10,10,10])*1e-3

pltext.initplot(num=1,xlabel='Frequenz in Hz',ylabel='Spannung in V',scale='loglog')
pltext.plotdata(f_1,U_A1,U_A1_dsys,label='R=680kOhm')
pltext.plotdata(f_1,U_A2,U_A2_dsys,label='R=274kOhm')
pltext.plotdata(f_1,U_A3,U_A3_dsys,label='R=48.7kOhm')
pltext.plotdata(f_1,U_A4,U_A4_dsys,label='mit C')

pltext.initplot(num=2,xlabel='Frequenz in Hz',ylabel='Spannung in V',scale='loglog')
pltext.plotdata(f_2,U_A5,U_A5_dsys)

plt.show()
