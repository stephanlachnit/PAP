# measure version 1.4.2
from measure import sqrt,log10,pi,mean_value,std_dev_m,plot,linreg,val,lst,sig

# measured values
s_stokes = [0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05]
t_stokes = [[8.90,8.98,8.90,8.95,8.93],[7.20,7.56,7.40,7.40,7.65],[8.26,8.20,8.18,8.19,7.98],[11.07,11.04,11.06,11.31,11.01],[15.39,15.18,14.61,14.62,14.96],[11.43,11.54,11.62,11.29,11.40],[20.07,19.40,19.75,19.39,19.36],[19.04,19.50,18.73,19.56,19.10],[31.54,33.43,33.18,32.93,32.43]]
V_hp = [0.0, 5e-6, 10e-6, 15e-6, 20e-6, 25e-6]
t_hp = [0.0, 2*60+14.21, 4*60+16.76, 6*60+37.78, 8*60+49.48, 11*60+05.76]
hA = 526e-3
dhA = 1e-3
hE = 520e-3
dhE = 1e-3
dt_systematic = sqrt(2) * 0.2
r = [9e-3, 8e-3, 7.144e-3, 6e-3, 5e-3, 4e-3, 3e-3, 2e-3, 1.5e-3]
dr = [0.01 * r[i] for i in range(len(r))]
d_s = [1362.5, 1357.5, 1377.5, 1377.5, 1377.5, 1377.5, 1377.5, 1377.5, 1392.5]
dd_s = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
d_f = 1147.6
dd_f = 1.1
R_stokes = 75e-3
dR_stokes = 1e-3
R_hp = 1/2 * 1.5e-3
dR_hp = 1/2 * 0.01e-3
L_hp = 100e-3
dL_hp = 0.5e-3
g = 9.80984
dg = 2e-5

# viscosity with stokes
N = len(t_stokes)
t = [mean_value(t_stokes[i]) for i in range(N)]
dt = [sqrt(std_dev_m(t_stokes[i])**2 + dt_systematic**2) for i in range(N)]

v_stokes = [s_stokes[i] / t[i] for i in range(len(s_stokes))]
dv_stokes = [s_stokes[i] / t[i]**2 * dt[i] for i in range(N)]

r2 = [r[i]**2 for i in range(N)]
dr2 = [2 * r[i] * dr[i] for i in range(N)]
vdd = [v_stokes[i] / (d_s[i] - d_f) for i in range(N)]
dvdd = [1 / (d_s[i] - d_f) * sqrt(dv_stokes[i]**2 + (v_stokes[i] / (d_s[i] - d_f))**2 * (dd_s[i]**2 + dd_f**2)) for i in range(N)]

corr = [1 + 2.1 * r[i] / R_stokes for i in range(N)]
dcorr = [2.1 / R_stokes * sqrt(dr[i]**2 + (r[i] * dR_stokes / R_stokes)**2) for i in range(N)]
vddcorr = [vdd[i] * corr[i] for i in range(N)]
dvddcorr = [sqrt((dvdd[i] * corr[i])**2 + (vdd[i] * dcorr[i])**2) for i in range(N)]

lrplot1 = plot(title="Comparison of the sphere's speeds with and without correction λ", xlabel="r^2 / m^2", ylabel="v/(ρs - ρf) / m^4*kg^-1*s^-1", figure=1)
[tmp, tmp, tmp, tmp] = linreg(r2, vddcorr, dvddcorr, dr2, drawplot=True, graphname="w/ λ", lrplot=lrplot1)
[tmp, tmp, tmp, tmp] = linreg(r2, vdd, dvdd, dr2, drawplot=True, graphname="w/o λ", lrplot=lrplot1)
lrplot1.drawplot()

r2_lr = [0.0, r2[8], r2[7], r2[6], r2[5], r2[4]]
dr2_lr = [0.0, dr2[8], dr2[7], dr2[6], dr2[5], dr2[4]]
vdd_lr = [0.0, vddcorr[8], vddcorr[7], vddcorr[6], vddcorr[5], vddcorr[4]]
dvdd_lr = [0.0, dvddcorr[8], dvddcorr[7], dvddcorr[6], dvddcorr[5], dvddcorr[4]]
lrplot2 = plot(title="Linear regression of the slope", xlabel="r^2 / m^2", ylabel="v/(ρs - ρf) / m^4*kg^-1*s^-1", figure=2)
[slope, dslope, tmp, tmp] = linreg(r2_lr, vdd_lr, dvdd_lr, dr2_lr, drawplot=True, lrplot=lrplot2)
lrplot2.drawplot()

n_stokes = 2 / 9 * g / slope
dn_stokes = 2 / 9 / slope * sqrt(dg**2 + (g / slope * dslope)**2)

Re_stokes = [d_f * v_stokes[i] * 2 * r[i] / n_stokes for i in range(N)]
dRe_stokes = [2 / n_stokes * sqrt((r[i] * v_stokes[i] * dd_f)**2 + (r[i] * d_f * dv_stokes[i])**2 + (d_f * v_stokes[i] * dr[i])**2 + (r[i] * d_f * v_stokes[i] / n_stokes * dn_stokes)**2) for i in range(N)]

vlam = [2 / 9 * g * (d_s[i] - d_f) / n_stokes * r2[i] for i in range(N)]
dvlam = [2 / 9 / n_stokes * sqrt((g * r2[i])**2 * (dd_s[i]**2 + dd_f**2 + ((dd_s[i] - dd_f) / n_stokes * dn_stokes)**2) + (r2[i] * (d_s[i] - d_f) * dg)**2 + (g * (dd_s[i] - dd_f) * dr2[i])**2) for i in range(N)]

vdvlam = [v_stokes[i] / vlam[i] for i in range(N)]
dvdlam = [1 / vlam[i] * sqrt(dv_stokes[i]**2 + (v_stokes[i] / vlam[i] * dvlam[i])**2) for i in range(N)]

lamplt = plot(xlabel="v / vlam", ylabel="Re", figure=3, scale="linlog")
lamplt.plotdata(vdvlam, Re_stokes, dRe_stokes, dvdlam)
lamplt.drawplot()

print()
print("viscosity with stokes:")
print()
print(lst("r", r, dr))
print()
print(lst("r2", r2, dr2))
print()
print(lst("v w/o λ", v_stokes, dv_stokes))
print()
print(lst("v/roh w/o λ", vdd, dvdd))
print()
print(lst("v/roh w/ λ", vddcorr, dvddcorr))
print()
print(val("slope", slope, dslope))
print()
print(val("viscosity", n_stokes, dn_stokes))
print()
print(lst("Re", Re_stokes, dRe_stokes))
print()
print(lst("vlam", vlam, dvlam))
print()
print(lst("v/vlam", vdvlam, dvdlam))

# viscosity with hagen-poiseuille
h = mean_value([hA, hE])
dh = sqrt(dhA**2 + dhE**2)
p = h * d_f * g
dp = sqrt((d_f * g * dh)**2 + (h * g * dd_f)**2 + (h * d_f *dg)**2)

hgplot = plot(title="Volume flow rate (Hagen-Poiseuille)", xlabel="V / m^3", ylabel="t / s", figure=4)
[tdV, dtdV, tmp, tmp] = linreg(V_hp, t_hp, [dt_systematic for i in range(len(t_hp))], drawplot=True, lrplot=hgplot)
hgplot.drawplot()

Vdt = 1 / tdV
dVdt = dtdV / tdV**2

n_hp = pi / 8 * p * R_hp**4 * tdV / L_hp
dn_hp = pi / 8 * R_hp**3 / L_hp * sqrt((p * 4 * dR_hp *tdV)**2 + (dp * R_hp * tdV)**2 + (p * R_hp * dtdV)**2 + (p * R_hp * tdV / L_hp * dL_hp)**2)

v_stokes = p * R_hp**2 / (6 * n_hp * h)
dv_hp = R_hp / (6 * n_hp * h) * sqrt((R_hp * dp)**2 + (R_hp * p * dn_hp / n_hp)**2 + (R_hp * p * dh / h)**2 + (2 * p * dR_hp)**2)
Re_hp = d_f * v_stokes * 2 * R_hp / n_hp
dRe_hp = 2 / n_hp * sqrt((R_hp * dd_f * v_stokes)**2 + (R_hp * d_f * dv_hp)**2 + (R_hp * d_f * v_stokes * dn_hp / n_hp)**2 + (dR_hp * d_f * v_stokes)**2)

print()
print("viscosity with hagen-poiseuille:")
print()
print(val("Volume flow speed", Vdt, dVdt))
print(val("Vertical flow speed", v_stokes, dv_hp))
print(val("Reynold", Re_hp, dRe_hp))
print(val("viscosity", n_hp, dn_hp))

# Comparison
print()
print(sig("Deviation n", n_hp, dn_hp, n_stokes, dn_stokes))

# Show Plots
plot.showplots()
