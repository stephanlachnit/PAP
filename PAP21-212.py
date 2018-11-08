# measure version 1.3
from measure import sqrt,T0,mean_value,std_dev_m,lst

# measured values
t_stokes = [[8.90,8.98,8.90,8.95,8.93],[7.20,7.56,7.40,7.40,7.65],[8.26,8.20,8.18,8.19,7.98],[11.07,11.04,11.06,11.31,11.01],[15.39,15.18,14.61,14.62,14.96],[11.43,11.54,11.62,11.29,11.40],[20.07,19.40,19.75,19.39,19.36],[19.04,19.50,18.73,19.56,19.10],[31.54,33.43,33.18,32.93,32.43]]
t_hp = [2*60+14.21, 4*60+16.76, 6*60+37.78, 8*60+49.48, 11*60+05.76]
dt_systematic = sqrt(2) * 0.2
s = [0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05]
r = [9e-3, 8e-3, 7.144e-3, 6e-3, 5e-3, 4e-3, 3e-3, 2e-3, 1.5e-3]
R = 75e-3
dR = 1e-3
T = 22.0 + T0
dT = 1.0

# viscosity with stokes
t = [mean_value(t_stokes[i]) for i in range(len(t_stokes))]
dt = [std_dev_m(t_stokes[i]) for i in range(len(t_stokes))]

v = [s[i] / t[i] for i in range(len(s))]
dv = [s[i] / t[i]**2 * dt[i] for i in range(len(s))]

print()
print("viscosity with stokes:")
print(lst("velocity", v, dv))
