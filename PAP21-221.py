# measure version 1.7.1
from measure import npfarray,sqrt,mv,dsto_mv,dsys_mv,val,sig,tbl

# values
cd_h1r = npfarray([0.610, 0.626, 0.569, 0.630, 0.570])
cd_h1l = npfarray([0.489, 0.471, 0.530, 0.467, 0.529])
cd_h3r = npfarray([0.565, 0.569, 0.553, 0.570, 0.556])
cd_h3l = npfarray([0.533, 0.530, 0.545, 0.529, 0.544])
cd_hix_dsys = npfarray([1e-3 for i in range(5)])

rh_50T_air = 50.04
rh_50T_arg = 46.57
rh_50T_dsys = 0.3
rh_V_air = 5370e-6
rh_V_arg = 5460e-6
rh_V_dsys = 5e-6
rh_m_air = 26.116e-3
rh_m_arg = 26.006e-3
rh_m_dsys = 0.002e-3
rh_2r_air = 15.95e-3
rh_2r_arg = 15.97e-3
rh_2r_dsys = 0.02e-3
rh_p = 994.15e2
rh_p_dsys = 0.15e2

k_air_lit = 0.78 * 1.401 + 0.21 * 1.398 + 0.01 * 1.648
k_arg_lit = 1.648

# Clément & Desormes
h1 = cd_h1r - cd_h1l
h3 = cd_h3r - cd_h3l
hi_dsys = sqrt(2) * cd_hix_dsys

cd_k = h1 / (h1 - h3)
cd_k_dsys = hi_dsys / (h1 - h3) * sqrt((1 + h1 / (h1 -h3))**2 + (1 / (h1 - h3))**2)

cd_k_mv = mv(cd_k)
cd_k_dsto_mv = dsto_mv(cd_k)
cd_k_dsys_mv = dsys_mv(cd_k_dsys)
cd_k_dtot = sqrt(cd_k_dsto_mv**2 + cd_k_dsys_mv**2)

print()
print('Clément & Desormes:')
print()
print(tbl(['h1','h3','k'], [h1, h3, cd_k], [hi_dsys, hi_dsys, cd_k_dsys]))
print()
print(val('k', cd_k_mv, cd_k_dtot))
print(sig('dev', cd_k_mv, cd_k_dtot, k_air_lit))

# Rüchardt
r_air = rh_2r_air / 2.
r_air_dsys = rh_2r_dsys / 2.
T_air = rh_50T_air / 50.
T_air_dsys = rh_50T_dsys / 50.
rh_k_air = 4. * rh_m_air * rh_V_air / (r_air**4 * T_air**2 * rh_p)
rh_k_air_dsys = 4. / (r_air**4 * T_air**2 * rh_p) * sqrt((rh_m_dsys * rh_V_air)**2 + (rh_m_air * rh_V_dsys)**2 + (rh_m_air * rh_V_air)**2 * ((4. * r_air_dsys / r_air)**2 + (2. * T_air_dsys / T_air)**2 + (rh_p_dsys / rh_p)**2))

r_arg = rh_2r_arg / 2.
r_arg_dsys = rh_2r_dsys / 2.
T_arg = rh_50T_arg / 50.
T_arg_dsys = rh_50T_dsys / 50.
rh_k_arg = (4. * rh_m_arg * rh_V_arg) / (r_arg**4 * T_arg**2 * rh_p)
rh_k_arg_dsys = 4. / (r_arg**4 * T_arg**2 * rh_p) * sqrt((rh_m_dsys * rh_V_arg)**2 + (rh_m_arg * rh_V_dsys)**2 + (rh_m_arg * rh_V_arg)**2 * ((4. * r_arg_dsys / r_arg)**2 + (2. * T_arg_dsys / T_arg)**2 + (rh_p_dsys / rh_p)**2))

print()
print('Rüchardt:')
print()
print('air:')
print(val('r',r_air, r_air_dsys))
print(val('T', T_air, T_air_dsys))
print(val('k', rh_k_air, rh_k_air_dsys))
print(sig('dev', rh_k_air, rh_k_air_dsys, k_air_lit))
print()
print('argon:')
print(val('r',r_arg, r_arg_dsys))
print(val('T', T_arg, T_arg_dsys))
print(val('k', rh_k_arg, rh_k_arg_dsys))
print(sig('dev', rh_k_arg, rh_k_arg_dsys, k_arg_lit))
