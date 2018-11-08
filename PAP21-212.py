# measure version 1.4
from measure import sqrt,log10,T0,mean_value,std_dev_m,plot,linreg,val,lst

# measured values
t_stokes = [[8.90,8.98,8.90,8.95,8.93],[7.20,7.56,7.40,7.40,7.65],[8.26,8.20,8.18,8.19,7.98],[11.07,11.04,11.06,11.31,11.01],[15.39,15.18,14.61,14.62,14.96],[11.43,11.54,11.62,11.29,11.40],[20.07,19.40,19.75,19.39,19.36],[19.04,19.50,18.73,19.56,19.10],[31.54,33.43,33.18,32.93,32.43]]
t_hp = [2*60+14.21, 4*60+16.76, 6*60+37.78, 8*60+49.48, 11*60+05.76]
dt_systematic = sqrt(2) * 0.2
s = [0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05]
r = [9e-3, 8e-3, 7.144e-3, 6e-3, 5e-3, 4e-3, 3e-3, 2e-3, 1.5e-3]
d_s = [1362.5, 1357.5, 1377.5, 1377.5, 1377.5, 1377.5, 1377.5, 1377.5, 1392.5]
dd_s = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
d_f = 1147.6
dd_f = 1.1
R = 75e-3
dR = 1e-3
T = 22.0 + T0
dT = 1.0
g = 9.80984
dg = 2e-5

# viscosity with stokes
N = len(t_stokes)
t = [mean_value(t_stokes[i]) for i in range(N)]
dt = [std_dev_m(t_stokes[i]) + dt_systematic for i in range(N)]

v = [s[i] / t[i] for i in range(len(s))]
dv = [s[i] / t[i]**2 * dt[i] for i in range(N)]

r2 = [r[i]**2 for i in range(N)]
vd = [v[i] / (d_s[i] - d_f) for i in range(N)]
dvd = [1 / (d_s[i] - d_f) * sqrt(dv[i]**2 + (v[i] / (d_s[i] - d_f))**2 * (dd_s[i]**2 + dd_f**2)) for i in range(N)]

corr = [1 + 2.1 * r[i] / R for i in range(N)]
vdcorr = [vd[i] * corr[i] for i in range(N)]
dvdcorr = [dvd[i] * corr[i] for i in range(N)]

lrplot1 = plot(title="all values", xlabel="r^2 / m^2", ylabel="v/(ρs - ρf) / m^2*kg^-1*s^-1", figure=1)
[tmp, tmp, tmp, tmp] = linreg(r2, vd, dvd, drawplot=True, graphname="w/o correction", lrplot=lrplot1)
[tmp, tmp, tmp, tmp] = linreg(r2, vdcorr, dvdcorr, drawplot=True, graphname="w/  correction", lrplot=lrplot1)

r2_lr = [0.0, r2[8], r2[7], r2[6], r2[5], r2[4], r2[3]]
vd_lr = [0.0, vdcorr[8], vdcorr[7], vdcorr[6], vdcorr[5], vdcorr[4], vdcorr[3]]
dvd_lr = [1e-6, dvdcorr[8], dvdcorr[7], dvdcorr[6], dvdcorr[5], dvdcorr[4], dvdcorr[3]]
lrplot2 = plot(title="used values", xlabel="r^2 / m^2", ylabel="v/(ρs - ρf) / m^2*kg^-1*s^-1", figure=2)
[slope, dslope, tmp, tmp] = linreg(r2_lr, vd_lr, dvd_lr, drawplot=True, lrplot=lrplot2)

n_stokes = 2 / 9 * g / slope
dn_stokes = 2 / 9 / slope * sqrt(dg**2 + (g / slope * dslope)**2)

Re = [d_f * v[i] * 2 * r[i] / n_stokes for i in range(N)]
dRe = [2 * r[i] / n_stokes * sqrt((v[i] * dd_f)**2 + (d_f * dv[i])**2 + (d_f * v[i] / n_stokes * dn_stokes)**2) for i in range(N)]

vlam = [2 / 9 * g * (d_s[i] - d_f) / n_stokes * r2[i] for i in range(N)]
dvlam = [2 / 9 * r2[i] / n_stokes * sqrt(g**2 * (dd_s[i]**2 + dd_f**2 + ((dd_s[i] - dd_f) / n_stokes * dn_stokes)**2) + ((d_s[i] - d_f) * dg)**2) for i in range(N)]

vdvlam = [v[i] / vlam[i] for i in range(N)]
dvdlam = [1 / vlam[i] * sqrt(dv[i]**2 + (v[i] / vlam[i] * dvlam[i]**2)) for i in range(N)]

lamplt = plot(xlabel="v / vlam", ylabel="Re", figure=3, scale="linlog")
lamplt.plotdata(vdvlam, Re, dRe, dvdlam)
lamplt.drawplot()

print()
print("viscosity with stokes:")
print(lst("velocity", v, dv))

plot.showplots()
