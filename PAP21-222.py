# measure version 1.8
from measure import sqrt,val,npfarray,mv,dsto_mv

# values
cw = 4180
rohw = 998 # https://www.internetchemie.info/chemie-lexikon/daten/w/wasser-dichtetabelle.php
rohw_dsys = 1

U_mot = 24.0
U_mot_dsys = 0.1
I_mot = 2.3
I_mot_dsys = 0.1

km_U_heiz = 5.58
km_U_heiz_dsys = 0.04
km_I_heiz = 5.7
km_I_heiz_dsys = 0.1
km_f = 325.2 / 60.
km_f_dsys = 0.1 / 60.
km_t = 180.
km_t_dsys = 15.
km_dT = 22.1 - 18.65
km_dT_dsys = sqrt(0.25**2 + 0.1**2)
Vpm = npfarray([199.5,198.4,200.7,200.1,200.5])
km_Vps = mv(Vpm / 6e5)
km_Vps_dsys = dsto_mv(Vpm / 6e5)

# Kältemaschine
Wm = U_mot * I_mot / km_f
Wm_dsys = 1./km_f * sqrt((U_mot * I_mot_dsys)**2 + (U_mot_dsys * I_mot)**2 + (U_mot * I_mot * km_f_dsys / km_f)**2)
Wh = km_U_heiz * km_I_heiz / km_f
Wh_dsys = 1./km_f * sqrt((km_U_heiz * km_I_heiz_dsys)**2 + (km_U_heiz_dsys * km_I_heiz)**2 + (km_U_heiz * km_I_heiz * km_f_dsys / km_f)**2)

km_n = Wh / Wm
km_n_dsys = 1./Wm * sqrt(Wh_dsys**2 + (Wh * Wm_dsys / Wm)**2)

tWges = cw * rohw * km_dT * km_Vps / km_f
tWges_dsys = cw / km_f * sqrt((rohw_dsys * km_dT * km_Vps)**2 + (rohw * km_dT_dsys * km_Vps)**2 + (rohw * km_dT * km_Vps_dsys)**2 + (rohw * km_dT * km_Vps * km_f_dsys / km_f)**2)
eWges = km_t  * km_f * (Wm + Wh)
eWges_dsys = sqrt((km_t_dsys  * km_f * (Wm + Wh))**2 + (km_t  * km_f_dsys * (Wm + Wh))**2 + (km_t  * km_f * (Wm_dsys + Wh_dsys))**2)

dW = eWges - tWges
dW_dsys = sqrt(tWges_dsys**2 + eWges_dsys**2)

print()
print(val('Wirkungsgrad', km_n, km_n_dsys))
print(val('Ernegiedifferenz', dW, dW_dsys))
