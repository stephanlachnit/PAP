# measure version 1.8.2
from measure import sqrt,pi,val,lst,tbl,npfarray,mv,dsto_mv,linreg,plt,pltext

# values
c_w = 4180.
lamba_w = 335e3
roh_w = 998.0
roh_w_dsys = 1.0

U_mot = 24.0
U_mot_dsys = 0.1
I_mot = 2.3
I_mot_dsys = 0.1

km_U_heiz = 5.58
km_U_heiz_dsys = 0.04
km_I_heiz = 1.14 * 5
km_I_heiz_dsys = 0.02 * 5
km_f = 325.2 / 60.
km_f_dsys = 0.1 / 60.
km_dT = 22.1 - 18.65
km_dT_dsys = sqrt(0.25**2 + 0.1**2)
km_Vps = npfarray([199.5,198.4,200.7,200.1,200.5]) / 6e7
km_Vps_mv = mv(km_Vps)
km_Vps_dsys = dsto_mv(km_Vps)

gf_t = 180.
gf_t_dsys = 15.
gf_V = 1e-6
gf_V_dsys = 0.5e-6
gf_f = 306.5 / 60.
gf_f_dsys = 0.1 / 60.

wk_l = 0.250
wk_l_dsys = 0.005
wk_Vps = 199.0 / 6e7
wk_Vps_dsys = 1.5 / 6e7
wk_dT = 26.3 - 19.0
wk_dT_dsys = sqrt(2) * 0.1
wk_F = npfarray([0.0,0.2,0.4,0.6,0.8])
wk_F_dsys = npfarray([0.0,0.01,0.01,0.01,0.01])
wk_Uh = npfarray([11.90,11.95,12.04,12.04,11.97])
wk_Uh_dsys = npfarray([0.05,0.05,0.02,0.05,0.05])
wk_Ih = npfarray([2.59,2.59,2.62,2.61,2.61]) * 5
wk_Ih_dsys = npfarray([0.02,0.02,0.02,0.02,0.02]) * 5
wk_f = [[344.2,344.1,344.4],[307.1,308.6,309.3],[288.8,289.1,291.0],[262.4,263.2,263.5],[233.1,231.2,237.6]]
wk_f_mv = npfarray([mv(wk_f[n])/60. for n in range(5)])
wk_f_dsto = npfarray([dsto_mv(wk_f[n])/60. for n in range(5)])
wk_Wpv = [[24295,24524,24390],[25623,26181,25791],[27035,27112,27100],[28287,27615,28339],[29393,29246,29123]]
wk_Wpv_mv = npfarray([mv(wk_Wpv[n])*1e-4 for n in range(5)])
wk_Wpv_dsto = npfarray([dsto_mv(wk_Wpv[n])*1e-4 for n in range(5)])

# Kältemaschine
Wm = U_mot * I_mot / km_f
Wm_dsys = 1./km_f * sqrt((U_mot * I_mot_dsys)**2 + (U_mot_dsys * I_mot)**2 + (U_mot * I_mot * km_f_dsys / km_f)**2)
Wh = km_U_heiz * km_I_heiz / km_f
Wh_dsys = 1./km_f * sqrt((km_U_heiz * km_I_heiz_dsys)**2 + (km_U_heiz_dsys * km_I_heiz)**2 + (km_U_heiz * km_I_heiz * km_f_dsys / km_f)**2)

tWges = c_w * roh_w * km_dT * km_Vps_mv / km_f
tWges_dsys = c_w / km_f * sqrt((roh_w_dsys * km_dT * km_Vps_mv)**2 + (roh_w * km_dT_dsys * km_Vps_mv)**2 + (roh_w * km_dT * km_Vps_dsys)**2 + (roh_w * km_dT * km_Vps_mv * km_f_dsys / km_f)**2)

eWges = Wm + Wh
eWges_dsys = sqrt(Wm_dsys**2 + Wh_dsys**2)
dW = eWges - tWges
dW_dsys = sqrt(eWges_dsys**2 + tWges_dsys**2)

km_n = Wh / Wm
km_n_dsys = 1./Wm * sqrt(Wh_dsys**2 + (Wh * Wm_dsys / Wm)**2)

print()
print(tbl([['Wm',val('',Wm,Wm_dsys)],['Wh',val('',Wh,Wh_dsys)],['Wm+Wh',val('',eWges,eWges_dsys)],['Q1',val('',tWges,tWges_dsys)],['dW',val('',dW,dW_dsys)]]))
print(val('Wirkungsgrad', km_n, km_n_dsys))

# Gefrierzeit Wasser
Pk = gf_V * roh_w * lamba_w / gf_t
Pk_dsys = lamba_w / gf_t * sqrt((gf_V * roh_w * gf_t_dsys / gf_t)**2 + (gf_V * roh_w_dsys)**2 + (gf_V_dsys * roh_w)**2)
Wk = Pk / gf_f
Wk_dsys = 1 / gf_f * sqrt(Pk_dsys**2 + (Pk * gf_f_dsys / gf_f)**2)

print()
print(val('Pk', Pk, Pk_dsys))
print(val('Wk', Wk, Wk_dsys))

# Wärmekraftmaschine
Qel = wk_Uh * wk_Ih / wk_f_mv
Qel_dtot = 1 / wk_f_mv * sqrt((wk_Uh * wk_Ih_dsys)**2 + (wk_Uh_dsys * wk_Ih)**2 + (wk_Uh * wk_Ih * wk_f_dsto / wk_f_mv)**2)
Qab = c_w * roh_w * wk_dT * wk_Vps / wk_f_mv
Qab_dtot = c_w / wk_f_mv * sqrt((roh_w_dsys * wk_dT * wk_Vps)**2 + (roh_w * wk_dT_dsys * wk_Vps)**2 + (roh_w * wk_dT * wk_Vps_dsys)**2 + (roh_w * wk_dT * wk_Vps * wk_f_dsto / wk_f_mv)**2)
W_D = 2.*pi * wk_l * wk_F
W_D_dsys = 2.*pi * sqrt((wk_l_dsys * wk_F)**2 + (wk_l * wk_F_dsys)**2)

Pel = Qel * wk_f_mv
Pel_dtot = sqrt((Qel_dtot * wk_f_mv)**2 + (Qab * wk_f_dsto)**2)
Pab = Qab * wk_f_mv
Pab_dtot = sqrt((Qab_dtot * wk_f_mv)**2 + (Qab * wk_f_dsto)**2)
Ppv = wk_Wpv_mv * wk_f_mv
Ppv_dsto = sqrt((wk_Wpv_dsto * wk_f_mv)**2 + (wk_Wpv_mv * wk_f_dsto)**2)
P_D = W_D * wk_f_mv
P_D_dtot = sqrt((W_D_dsys * wk_f_mv)**2 + (W_D * wk_f_dsto)**2)

Q_V_pv = Qel - Qab - wk_Wpv_mv
Q_V_pv_dtot = sqrt(Qel_dtot**2 + Qab_dtot**2 + wk_Wpv_dsto**2)
Q_V_D = Qel - Qab - W_D
Q_V_D_dtot = sqrt(Qel_dtot**2 + Qab_dtot**2 + W_D_dsys**2)

P_V_pv = Q_V_pv * wk_f_mv
P_V_pv_dtot = sqrt((Q_V_pv_dtot * wk_f_mv)**2 + (Q_V_pv * wk_f_dsto)**2)
P_V_D = Q_V_D * wk_f_mv
P_V_D_dtot = sqrt((Q_V_D_dtot * wk_f_mv)**2 + (Q_V_D * wk_f_dsto)**2)

wk_n_th = wk_Wpv_mv / Qel
wk_n_th_dtot = 1 / Qel * sqrt(wk_Wpv_dsto**2 + (wk_Wpv_mv * Qel_dtot / Qel)**2)
wk_n_eff = W_D / Qel
wk_n_eff_dsys = 1 / Qel * sqrt(W_D_dsys**2 + (W_D * Qel_dtot / Qel)**2)

print()
print(tbl([['Qel']+lst(Qel,Qel_dtot),['Qab']+lst(Qab,Qab_dtot),['Wpv']+lst(wk_Wpv_mv,wk_Wpv_dsto),['W_D']+lst(W_D,W_D_dsys)]))
print(tbl([['Pel']+lst(Pel,Pel_dtot),['Pab']+lst(Pab,Pab_dtot),['Ppv']+lst(Ppv,Ppv_dsto),['P_D']+lst(P_D,P_D_dtot)]))
print(tbl([['Q_V (Wpv)']+lst(Q_V_pv,Q_V_pv_dtot),['Q_V (W_D)']+lst(Q_V_D,Q_V_D_dtot),['P_V (Wpv)']+lst(P_V_pv,P_V_pv_dtot),['P_V (W_D)']+lst(P_V_D,P_V_D_dtot)]))
print(tbl([['f']+lst(wk_f_mv,wk_f_dsto),['F']+lst(wk_F,wk_F_dsys),['n_th']+lst(wk_n_th,wk_n_th_dtot),['n_eff']+lst(wk_n_eff,wk_n_eff_dsys)]))

pltext.initplot(title='Wirkungsgrade in Abhängigkeit von der Frequenz', xlabel='Frequenz f / Hz', ylabel='Wirkungsgrad')
pltext.plotdata(wk_f_mv, wk_n_th, wk_n_th_dtot, wk_f_dsto, label=r'$n_{th}$', connect=True)
pltext.plotdata(wk_f_mv, wk_n_eff, wk_n_eff_dsys, wk_f_dsto, label=r'$n_{eff}$', connect=True)
plt.xlim(3.75,5.875)
plt.ylim(0.0,0.1)
plt.legend(loc='upper left')
plt.savefig('fig0.pdf', format='pdf')
plt.show()
