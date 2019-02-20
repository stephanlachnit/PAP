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

pltext.initplot(num=1,title='Phase',scale='loglin')
pltext.plotdata(f,Phi,Phi_dsys)


# Aufgabe 5
A = npfarray([1.58,1.13,0.80,0.58,0.41])
A_dsys = npfarray([0.05,0.05,0.05,0.05,0.05])
T = npfarray([260,260,257,257,258])
T_dsys = npfarray([10,10,5,5,5])
T_mv = mv(T)
T_mv_dtot = dtot(dsys_mv(T_dsys),dsto_mv(T))

t = npfarray([n*T_mv for n in range(0,5)])
dt = npfarray([T_mv_dtot for n in range(0,5)])

pltext.initplot(num=2,title='Dämpfung',scale='linlog')
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
plt.show()
