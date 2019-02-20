# measure version: 1.7.1
from measure import np,pi,g,dg,ln,exp,sqrt,val,lst,tbl,sig,mv,dsto_mv,dsys_mv,dtot,linreg,plot

# values
m_k = 4.164
m_s = 9.85e-3
r_k = 50.8e-3

damping_umin = np.array([1070, 970, 890, 808, 737, 675, 620, 569], dtype='float')
damping_dumin = np.array([10, 5, 4, 2, 3, 3, 2, 4], dtype='float')
damping_t = np.array([0, 2*60, 4*60, 6*60, 8*60, 10*60, 12*60, 14*60], dtype='float')

moa_umin_1x15 = np.array([640, 560, 420, 300], dtype='float')
moa_Tp_1x15 = np.array([104, 90, 70, 49], dtype='float')
moa_umin_1x20 = np.array([700, 600, 430, 270], dtype='float')
moa_Tp_1x20 = np.array([89, 77, 57, 36], dtype='float')
moa_umin_2x15 = np.array([670, 580, 430, 300], dtype='float')
moa_Tp_2x15 = np.array([60, 52, 39, 28], dtype='float')
moa_umin_2x20 = np.array([680, 540, 430, 320], dtype='float')
moa_Tp_2x20 = np.array([47, 38, 29, 21], dtype='float')
moa_dumin = np.array([10 for i in range(4)], dtype='float')
moa_dTp = np.array([1 for i in range(4)], dtype='float')

cf_f = np.array([570, 610, 470, 425, 400, 370, 320, 490, 540, 650], dtype='float')
cf_df = np.array([5 for i in range(10)], dtype='float')
cf_10u = np.array([17, 16, 20, 22, 23, 25, 29, 20, 17, 15], dtype='float')
cf_d10u = np.array([1 for i in range(10)], dtype='float')

nut_umin_n = np.array([670, 1025, 570, 300, 420, 475, 655, 720, 380, 330], dtype='float')
nut_dumin_n = np.array([10, 25, 10, 10, 10, 10, 10, 10, 10, 10], dtype='float')
nut_umin_f = np.array([720, 1120, 620, 320, 450, 520, 710, 770, 410, 360], dtype='float')
nut_dumin_f = np.array([10, 25, 10, 10, 10, 10, 10, 10, 10, 10], dtype='float')

# damping
f = damping_umin / 60.
df = damping_dumin / 60.
t = damping_t

[slope, dslope, yitc, dyitc] = linreg(x=t, y=ln(f), dy=df/f)
D = -1. * slope
dD = dslope
f0 = exp(yitc)
df0 = exp(yitc) * dyitc

t_hv = ln(2.) / D
t_hv_dsys = ln(2.) * dD / D**2

dampingplot = plot(title='Linear regression of the damping constant', xlabel='t / s', ylabel='f / Hz', fig=1, scale='linlog')
dampingplot.plt.grid(False)
dampingplot.plotdata(t, f, df)
dampingplot.plotfunc(t, f0 * exp(-D * t), label='line of best fit')
dampingplot.plotfunc(t, (f0+df0) * exp(-(D+dD) * t), label='line of uncertanty')

print()
print(val('damping constant', D, dD))
print(val('half value time', t_hv, t_hv_dsys))

# moment of anertia Iz with precession
f_m1x15 = moa_umin_1x15 / 120. * (1. + exp(-D * moa_Tp_1x15))
df_m1x15 = 1. / 120. * sqrt((moa_dumin * (1. + exp(-D * moa_Tp_1x15)))**2 + (moa_umin_1x15 * exp(-D * moa_Tp_1x15))**2 * ((D * moa_dTp)**2 + (dD * moa_Tp_1x15)**2))
f_m1x20 = moa_umin_1x20 / 120. * (1. + exp(-D * moa_Tp_1x20))
df_m1x20 = 1. / 120. * sqrt((moa_dumin * (1. + exp(-D * moa_Tp_1x20)))**2 + (moa_umin_1x20 * exp(-D * moa_Tp_1x20))**2 * ((D * moa_dTp)**2 + (dD * moa_Tp_1x20)**2))
f_m2x15 = moa_umin_2x15 / 120. * (1. + exp(-D * moa_Tp_2x15))
df_m2x15 = 1. / 120. * sqrt((moa_dumin * (1. + exp(-D * moa_Tp_2x15)))**2 + (moa_umin_2x15 * exp(-D * moa_Tp_2x15))**2 * ((D * moa_dTp)**2 + (dD * moa_Tp_2x15)**2))
f_m2x20 = moa_umin_2x20 / 120. * (1. + exp(-D * moa_Tp_2x20))
df_m2x20 = 1. / 120. * sqrt((moa_dumin * (1. + exp(-D * moa_Tp_2x20)))**2 + (moa_umin_2x20 * exp(-D * moa_Tp_2x20))**2 * ((D * moa_dTp)**2 + (dD * moa_Tp_2x20)**2))

moaplot = plot(title='precession period as function of the rotation frequency', xlabel='f / Hz', ylabel=r'$T_p$ / s', fig=2)
[slope_1x15, dslope_1x15, yitc_1x15, dyitc_1x15] = linreg(f_m1x15, moa_Tp_1x15, moa_dTp, df_m1x15, lrplot=moaplot, graphname='1@15')
[slope_1x20, dslope_1x20, yitc_1x20, dyitc_1x20] = linreg(f_m1x20, moa_Tp_1x20, moa_dTp, df_m1x20, lrplot=moaplot, graphname='1@20')
[slope_2x15, dslope_2x15, yitc_2x15, dyitc_2x15] = linreg(f_m2x15, moa_Tp_2x15, moa_dTp, df_m2x15, lrplot=moaplot, graphname='2@15')
[slope_2x20, dslope_2x20, yitc_2x20, dyitc_2x20] = linreg(f_m2x20, moa_Tp_2x20, moa_dTp, df_m2x20, lrplot=moaplot, graphname='2@20')

Iz_1x15 = m_s * g * 0.15 * slope_1x15 / (2.*pi)**2
dIz_1x15 = m_s * 0.15 / (2.*pi)**2 * sqrt((g * dslope_1x15)**2 + (dg * slope_1x15)**2)
Iz_1x20 = m_s * g * 0.20 * slope_1x20 / (2.*pi)**2
dIz_1x20 = m_s * 0.20 / (2.*pi)**2 * sqrt((g * dslope_1x20)**2 + (dg * slope_1x20)**2)
Iz_2x15 = 2.*m_s * g * 0.15 * slope_2x15 / (2.*pi)**2
dIz_2x15 = 2.*m_s * 0.15 / (2.*pi)**2 * sqrt((g * dslope_2x15)**2 + (dg * slope_2x15)**2)
Iz_2x20 = 2.*m_s * g * 0.20 * slope_2x20 / (2.*pi)**2
dIz_2x20 = 2.*m_s * 0.20 / (2.*pi)**2 * sqrt((g * dslope_2x20)**2 + (dg * slope_2x20)**2)

Iz_list = np.array([Iz_1x15, Iz_1x20, Iz_2x15, Iz_2x20])
dIz_list = np.array([dIz_1x15, dIz_1x20, dIz_2x15, dIz_2x20])
Iz = mv(Iz_list)
Iz_dsto = dsto_mv(Iz_list)
Iz_dsys = dsys_mv(dIz_list)
Iz_dtot = dtot(Iz_dsto, Iz_dsys)

tblstr = ['1@15','1@20','2@15','2@20']

print()
print('Frequency f:')
print(tbl(tblstr, [f_m1x15, f_m1x20, f_m2x15, f_m2x20], [df_m1x15, df_m1x20, df_m2x15, df_m2x20]))
print()
print('Linreg results (slope / yitc):')
print(tbl(tblstr, [[slope_1x15, yitc_1x15], [slope_1x20, yitc_1x20], [slope_2x15, yitc_2x15], [slope_2x20, yitc_2x20]], [[dslope_1x15, dyitc_1x15], [dslope_1x20, dyitc_1x20], [dslope_2x15, dyitc_2x15], [dslope_2x20, dyitc_2x20]]))
print()
print(sig('0-yitc 1@15', yitc_1x15, dyitc_1x15, 0.0))
print(sig('0-yitc 1@20', yitc_1x20, dyitc_1x20, 0.0))
print(sig('0-yitc 2@15', yitc_2x15, dyitc_2x15, 0.0))
print(sig('0-yitc 2@20', yitc_2x20, dyitc_2x20, 0.0))
print()
print(lst('Moment of anertia Iz', Iz_list, dIz_list))
print()
print(val('Iz (mean value)', Iz, Iz_dtot))

# color frequency
Omega = 2.*pi * 10. / cf_10u
dOmega = 2.*pi * 10. * cf_d10u / cf_10u**2
wf = 2.*pi * cf_f
dwf = 2.*pi * cf_df

cfplot = plot(title=r'angular color speed Ω as function of the angular rotation speed $ω_f$', xlabel=r'$ω_f$ / rad/s', ylabel='Ω / rad/s', fig=3)
[slope, dslope, yitc, dyitc] = linreg(wf, Omega, dOmega, dwf, lrplot=cfplot)

cf_Ix = Iz / (1 - slope)
cf_dIx = 1 / (1 - slope) * sqrt(Iz_dtot**2 + (Iz * dslope / (1 - slope))**2)

print(val('Ix (color)     ', cf_Ix, cf_dIx))

# nutation
f_f = nut_umin_f / 60.
df_f = nut_dumin_f / 60.
f_n = nut_umin_n / 60.
df_n = nut_dumin_n / 60.

nutplot = plot(title = r'nutation frequency $f_n$ as a function of figure frequency $f_f$', xlabel=r'$f_f$ / Hz', ylabel=r'$f_n$ / Hz', fig=4)
[slope, dslope, yitc, dyitc] = linreg(f_f, f_n, df_n, df_f, lrplot=nutplot)

nut_Ix = Iz / slope
nut_dIx = 1 / slope * sqrt(Iz_dtot**2 + (Iz * dslope / slope)**2)

print(val('Ix (nutation)  ', nut_Ix, nut_dIx))

# comparison
I_sphere = 2./5. * m_k * r_k**2

print()
print(sig('deviation Ix', cf_Ix, cf_dIx, nut_Ix, nut_dIx))
print()
print(val('Theoretical Value', I_sphere))
print(sig('dev Iz', Iz, Iz_dtot, I_sphere))
print(sig('dev Ix (cf)', cf_Ix, cf_dIx, I_sphere))
print(sig('dev Ix (nut)', nut_Ix, nut_dIx, I_sphere))
print(val('Iz / Ix (cf)', Iz / cf_Ix))
print(val('Iz / Ix (nut)', Iz / nut_Ix))

# show matplotlib figs
print()
plot.showfigs()
