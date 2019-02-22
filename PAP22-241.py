# measure version 1.8.7
from measure import npfarray,sqrt,ln,exp,arctan,lst,tbl,sig,val,mv,dsto_mv,dsys_mv,dtot_mv,plt,pltext,expreg,pi,curve_fit

# Aufgabe 1
R_A1 = npfarray([1,10,1])*1e3
R_A1_dsys = 0.05 * R_A1
C_A1 = npfarray([470,4.7,47])*1e-9
C_A1_dsys = 0.10 * C_A1
g_thalb = npfarray([312,32.6,32.6])*1e-6
g_thalb_dsys = npfarray([4,0.6,0.6])*1e-6

tau = R_A1 * C_A1
tau_dsys = sqrt((R_A1 * C_A1_dsys)**2 + (R_A1_dsys * C_A1)**2)
b_thalb = ln(2) * tau
b_thalb_dsys = ln(2) * tau_dsys

print()
print('Aufgabe 1:\n')
print(tbl([lst(R_A1,R_A1_dsys,'R'),lst(C_A1,C_A1_dsys,'C'),lst(tau,tau_dsys,'Tau')]))
print(tbl([lst(b_thalb,b_thalb_dsys,'T_1/2 (b)'),lst(g_thalb,g_thalb_dsys,'T_1/2 (g)'),['Abw']+[sig('',b_thalb[i],b_thalb_dsys[i],g_thalb[i],g_thalb_dsys[i]) for i in range(len(b_thalb))]]))

# Aufgabe 3
R_A3 = 1e3
R_A3_dsys = 0.05 * R_A3
C_A3 = 47e-9
C_A3_dsys = 0.05 * C_A3

f_A3 = npfarray([1,2,3,3.58,4,5,6,7,8,9,10])*1e3
dt = npfarray([200,83,46,30,27.2,18.8,14.0,10.4,8.0,6.8,4.1])*1e-6
dt_dsys = npfarray([20,5,5,5,3.0,2.5,2.5,2.0,1.5,1.5,1.0])*1e-6
Phi = 2*pi * f_A3 * dt
Phi_dsys = 2*pi * f_A3 * dt_dsys

from numpy import linspace
from scipy.optimize import fsolve

def phase_b(f):
  return arctan(1/(2*pi*f*R_A3*C_A3))
def phase_b_dys(f):
  return arctan(1/(2*pi*f*(R_A3+R_A3_dsys)*(C_A3+C_A3_dsys)))
def phase_b_45deg(f):
  return phase_b(f) - pi/4
def phase_b_45deg_dsys(f):
  return phase_b_dys(f) - pi/4

fgr_phase = fsolve(phase_b_45deg,x0=3.4e3)[0]
fgr_phase_dsys = abs(fsolve(phase_b_45deg_dsys,x0=0.3e3)[0] - fgr_phase)
fgr_fgang = npfarray([3.16,3.58])*1e3
fgr_fgang_dsys = npfarray([0.15,0.15])*1e3
fgr_fgang_mv = mv(fgr_fgang)
fgr_fgang_mv_dtot = dtot_mv(fgr_fgang,fgr_fgang_dsys)
fgr_calc = 1/(2*pi * R_A3 * C_A3)
fgr_calc_dsys = 1/(2*pi * R_A3 * C_A3) * sqrt((R_A3_dsys/R_A3)**2 + (C_A3_dsys/C_A3)**2)

f_array = linspace(1e3,10e3,1000)

pltext.initplot(num=1,title='Abbildung   : Phase in Abhängigkeit der Frequenz',xlabel='Frequenz in Hz',ylabel='Phase in rad',scale='loglin')
pltext.plotdata(f_A3,Phi,Phi_dsys,label='gemessene Phase')
plt.plot([1e3,10e3],[pi/4,pi/4],label='45°')
plt.plot(f_array,phase_b(f_array),label='berechnet')
plt.plot(f_array,phase_b_dys(f_array),label='berechnet, Fehler')
plt.legend()

print()
print('Aufgabe 3:\n')
print(tbl([['Messgröße','bei 45° Phase','Frequenzgang','berechnet'],lst([fgr_phase,fgr_fgang_mv,fgr_calc],[fgr_phase_dsys,fgr_fgang_mv_dtot,fgr_calc_dsys],'f_gr in Hz')]))
print(sig('Phase/Fgang',fgr_phase,fgr_phase_dsys,fgr_fgang_mv,fgr_fgang_mv_dtot))
print(sig('Phase/calc ',fgr_phase,fgr_phase_dsys,fgr_calc,fgr_calc_dsys))
print(sig('Fgang/calc ',fgr_fgang_mv,fgr_fgang_mv_dtot,fgr_calc,fgr_calc_dsys))

# Aufgabe 4
f_R = npfarray([3.93,3.73,3.70])*1e3
f_R_dsys = npfarray([0.10,0.05,0.03])*1e3
df = npfarray([4.90,1.26,0.57])*1e3
df_dsys = npfarray([0.05,0.05,0.05])*1e3
Ue = npfarray([0.99,0.96,0.95])
Ua = npfarray([0.95,0.77,0.35])
Ua_dsys = npfarray([0.01,0.01,0.01])

C_A4 = 47e-9
C_A4_dsys = 0.10 * C_A4
R_A4 = npfarray([1e3,220,47])
R_A4_dsys = 0.05 * R_A4

L_A4 = 1/(C_A4 * (2*pi * f_R)**2)
L_A4_dsys = 1/(C_A4 * (2*pi * f_R)**2) * sqrt((C_A4_dsys/C_A4)**2 + (2*f_R_dsys/f_R)**2)
L_A4_mv = mv(L_A4)
L_A4_mv_dtot = dtot_mv(L_A4,L_A4_dsys)

R_ges_A4 = 2*pi * df * L_A4
R_ges_A4_dsys = 2*pi * sqrt((df * L_A4_dsys)**2 + (df_dsys * L_A4)**2) 
R_ges_A4_mv = mv(R_ges_A4)
R_ges_A4_mv_dtot = dtot_mv(R_ges_A4,R_ges_A4_dsys)

R_V_df_A4 = R_ges_A4 - R_A4
R_V_df_A4_dsys = sqrt(R_ges_A4_dsys**2 + R_A4_dsys**2)

R_V_df_A4_mv = mv(R_V_df_A4)
R_V_df_A4_mv_dtot = dtot_mv(R_V_df_A4,R_ges_A4_dsys)

R_V_U_A4 = R_A4 * (Ue / Ua - 1)
R_V_U_A4_dsys = sqrt(((Ue / Ua - 1) * R_A4_dsys)**2 + (R_A4 * Ue / Ua**2 * Ua_dsys)**2)

R_V_U_A4_mv = mv(R_V_U_A4)
R_V_U_A4_mv_dtot = dtot_mv(R_V_U_A4,R_V_U_A4_dsys)

print()
print('Aufgabe 4:\n')
print(val(L_A4_mv,L_A4_mv_dtot,'Induktivität L'))
print()
print(tbl([lst(R_A4,R_A4_dsys,'R'),lst(R_ges_A4,R_ges_A4_dsys,'R_ges (df)'),lst(R_V_df_A4,R_V_df_A4_dsys,'R_V (df)'),lst(R_V_U_A4,R_V_U_A4_dsys,'R_V (U)')]))
print(tbl([lst([R_V_df_A4_mv],[R_V_df_A4_mv_dtot],'mv(R_V) (df)'),lst([R_V_U_A4_mv],[R_V_U_A4_mv_dtot],'mv(R_V) (U)'),['Abw']+[sig('',R_V_df_A4_mv,R_V_df_A4_mv_dtot,R_V_U_A4_mv,R_V_U_A4_mv_dtot)]]))

# Aufgabe 5
A = npfarray([1.58,1.13,0.80,0.58,0.41])
A_dsys = npfarray([0.05,0.05,0.05,0.05,0.05])
T = npfarray([260,260,257,257,258])*1e-6
T_dsys = npfarray([10,10,5,5,5])*1e-6
T_mv = mv(T)
T_mv_dtot = dtot_mv(T,T_dsys)

t = npfarray([n*T_mv for n in range(0,5)])
dt = npfarray([T_mv_dtot for n in range(0,5)])

pltext.initplot(num=2,title='Abbildung   : Bestimmung der Dämpfungskonstante',xlabel='Zeit in s',ylabel='Amplitude in V',scale='linlog')
dc,dc_dsys,yitc,dyitc = expreg(t,A,A_dsys,dt,plot=True)
dc *= -1

R_A5 = 47.
R_A5_dsys = 0.05 * R_A5
R_ges_A5 = dc * 2 * L_A4_mv
R_ges_A5_dsys = 2* sqrt((dc * L_A4_mv_dtot)**2 + (dc_dsys * L_A4_mv)**2)
R_V_A5 = R_ges_A5 - R_A5
R_V_A5_dsys = sqrt(R_ges_A5_dsys**2 + R_A5_dsys**2)

print()
print('Aufgabe 5:\n')
print(val(dc,dc_dsys,'Dämpfungskonstante d'))
print(val(R_ges_A5,R_ges_A5_dsys,'Gesamtwiderstand R'))
print()
print(tbl([lst([R_V_A5],[R_V_A5_dsys],'R_V (A5)'),['Abw df',sig('',R_V_A5,R_V_A5_dsys,R_V_df_A4_mv,R_V_df_A4_mv_dtot)],['Abw U',sig('',R_V_A5,R_V_A5_dsys,R_V_U_A4_mv,R_V_U_A4_mv_dtot)]]))

# Aufgabe 6
R_A6 = 220.
R_A6_dsys = 0.05 * R_A6
C_A6 = 47e-9
C_A6_dsys = 0.05 * C_A6
f_R_g = npfarray([3.75,3.94,3.85])*1e3
f_R_g_dsys = npfarray([0.03,0.03,0.05])*1e3

f_R_b = npfarray([0,0,0])
f_R_b_dtot = npfarray([0,0,0])
wr = 1/sqrt(L_A4_mv * C_A6)
wr_dtot = 1/sqrt(L_A4_mv * C_A6) * sqrt((C_A6_dsys / (2 * C_A6))**2 + (L_A4_mv_dtot / (2 * L_A4_mv))**2)
delta = R_A6 / (2 * L_A4_mv)
delta_dtot = 1/(2 * L_A4_mv) * sqrt(R_A6_dsys**2 + (R_A6 * L_A4_mv_dtot / L_A4_mv)**2)
f_R_b[0] = sqrt(wr**2 - 2 * delta**2)
f_R_b_dtot[0] = 1/f_R_b[0] * sqrt((wr * wr_dtot)**2 + (2 * delta * delta_dtot)**2)
f_R_b[1] = sqrt(wr**2 + 2 * delta**2)
f_R_b_dtot[1] = 1/f_R_b[1] * sqrt((wr * wr_dtot)**2 + (2 * delta * delta_dtot)**2)
f_R_b[2] = wr
f_R_b_dtot[2] = wr_dtot
f_R_b /= 2*pi
f_R_b_dtot /= 2*pi

print()
print('Aufgabe 6:\n')
print(tbl([ ['Größe','Kondensator','Spule','Widerstand'], lst(f_R_g,f_R_g_dsys,'f_R (g)'), lst(f_R_b,f_R_b_dtot,'f_R (b)'), ['Abw']+[sig('',f_R_g[i],f_R_g_dsys[i],f_R_b[i],f_R_b_dtot[i]) for i in range(0,3)] ]))

# Aufgabe 7
f_R_g = 3.89e3
f_R_g_dsys = 0.05e3

print()
print('Aufgabe 7:\n')
print(tbl([ lst([f_R_g],[f_R_g_dsys],'f_R (g)'), lst([f_R_b[2]],[f_R_b_dtot[2]],'f_R (b)'), ['Abw',sig('',f_R_g,f_R_g_dsys,f_R_b[2],f_R_b_dtot[2])] ]))

# Aufgabe 8
def V(dbV):
  return 10**(dbV/20)

dbV_vals = npfarray([[-2.8,-10.9,-20.0],[-32.3,-13.6,-20.8],[-2.8,-14.7,-27.8],[-2.7,8.3,-30.5],[-32.2,-11.3,-23.8]])
V_vals = V(dbV_vals)

perc_100 = [V_vals[i][0] / V_vals[0][0] for i in range(1,5)]
perc_4k = [V_vals[i][1] / V_vals[0][1] for i in range(1,5)]
perc_8k = [V_vals[i][2] / V_vals[0][2] for i in range(1,5)]

print()
print('Aufgabe 8:\n')
print('Spannungen in Volt')
print(tbl([ ['Signal','ohne Filter','RC-Hochpass','RC-Tiefpass','LC-Tiefpass','Bandpass 1k'], lst(V_vals[:,0],name='100Hz'), lst(V_vals[:,1],name='4kHz'), lst(V_vals[:,2],name='8kHz') ]))
print()
print('Verhältnis zum ungefilterten Signal')
print(tbl([ ['Signal','RC-Hochpass','RC-Tiefpass','LC-Tiefpass','Bandpass 1k'], lst(perc_100,name='100Hz'), lst(perc_4k,name='4kHz'), lst(perc_8k,name='8kHz') ]))

# Plot
print()
plt.show()
