# measure version 1.8.7
from measure import npfarray,sqrt,ln,lst,tbl,sig,val,mv,dsto_mv,dsys_mv,dtot,plt,pltext,expreg,pi

# Aufgabe 1
R = npfarray([1,10,1])*1e3
R_dsys = 0.05 * R
C = npfarray([470,4.7,47])*1e-9
C_dsys = 0.10 * C
g_thalb = npfarray([312,32.6,32.6])*1e-6
g_thalb_dsys = npfarray([4,0.6,0.6])*1e-6

tau = R * C
tau_dsys = sqrt((R * C_dsys)**2 + (R_dsys * C)**2)
b_thalb = ln(2) * tau
b_thalb_dsys = ln(2) * tau_dsys

print()
print('Aufgabe 1:')
print(tbl([lst(R,R_dsys,'R'),lst(C,C_dsys,'C'),lst(tau,tau_dsys,'Tau')]))
print(tbl([lst(b_thalb,b_thalb_dsys,'T_1/2 (b)'),lst(g_thalb,g_thalb_dsys,'T_1/2 (g)'),['Abw']+[sig('',b_thalb[i],b_thalb_dsys[i],g_thalb[i],g_thalb_dsys[i]) for i in range(len(b_thalb))]]))

# Aufgabe 3
f = npfarray([1,2,3,3.58,4,5,6,7,8,9,10])*1e3
dt = npfarray([200,83,46,30,27.2,18.8,14.0,10.4,8.0,6.8,4.1])*1e-6
dt_dsys = npfarray([20,5,5,5,3.0,2.5,2.5,2.0,1.5,1.5,1.0])*1e-6
Phi = 2*pi * f * dt
Phi_dsys = 2*pi * f * dt_dsys

pltext.initplot(num=1,title='Abbildung   : Phase in Abhängigkeit der Frequenz',xlabel='Frequenz in Hz',ylabel='Phase in rad',scale='loglin')
pltext.plotdata(f,Phi,Phi_dsys,label='gemessene Phase')
plt.plot([1e3,10e3],[pi/4,pi/4],label='45°')
plt.legend()

R_A3 = 1e3
R_A3_dsys = 0.05 * R_A3
C_A3 = 47e-9
C_A3_dsys = 0.05 * C_A3

fgr_phase = 3400
fgr_phase_dsys = 300
fgr_fgang = npfarray([3.16,3.58])*1e3
fgr_fgang_dsys = npfarray([0.15,0.15])*1e3
fgr_fgang_mv = mv(fgr_fgang)
fgr_fgang_dtot = dtot(dsys_mv(fgr_fgang_dsys),dsto_mv(fgr_fgang))
fgr_calc = 1/(2*pi * R_A3 * C_A3)
fgr_calc_dsys = 1/(2*pi * R_A3 * C_A3) * sqrt((R_A3_dsys/R_A3)**2 + (C_A3_dsys/C_A3)**2)

print()
print('Aufgabe 3:')
print(tbl([['','bei 45° Phase','Frequenzgang','berechnet'],lst([fgr_phase,fgr_fgang_mv,fgr_calc],[fgr_phase_dsys,fgr_fgang_dtot,fgr_calc_dsys],'f_gr in Hz')]))
print(sig('Phase/Fgang',fgr_phase,fgr_phase_dsys,fgr_fgang_mv,fgr_fgang_dtot))
print(sig('Phase/calc ',fgr_phase,fgr_phase_dsys,fgr_calc,fgr_calc_dsys))
print(sig('Fgang/calc ',fgr_fgang_mv,fgr_fgang_dtot,fgr_calc,fgr_calc_dsys))

# Aufgabe 4


# Aufgabe 5
A = npfarray([1.58,1.13,0.80,0.58,0.41])
A_dsys = npfarray([0.05,0.05,0.05,0.05,0.05])
T = npfarray([260,260,257,257,258])
T_dsys = npfarray([10,10,5,5,5])
T_mv = mv(T)
T_mv_dtot = dtot(dsys_mv(T_dsys),dsto_mv(T))

t = npfarray([n*T_mv for n in range(0,5)])
dt = npfarray([T_mv_dtot for n in range(0,5)])

pltext.initplot(num=2,title='Abbildung   : Bestimmung der Dämpfungskonstante',scale='linlog')
dc,dc_dsys,yitc,dyitc = expreg(t,A,A_dsys,dt,plot=True)

L = L_dsys = 1
R_ges = dc * 2 * L
R_ges_dsys = 2* sqrt((dc * L_dsys)**2 + (dc_dsys * L)**2)

print()
print('Aufgabe 5:')
print(val(dc,dc_dsys,'Dämpfungskonstante'))
print(val(R_ges,R_ges_dsys,'R'))

# Plot
print()
#plt.show()
