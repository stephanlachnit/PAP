# measure version: 1.5.1
from measure import np,ln,exp,linreg,plot,showfigs

# values
damping_umin = np.array([1070, 970, 890, 808, 737, 675, 620, 569], dtype='float')
damping_dumin = np.array([10, 5, 4, 2, 3, 3, 2, 4], dtype='float')
damping_t = np.array([0, 2, 4, 6, 8, 10, 12, 14], dtype='float')

# damping
f = damping_umin * 60
df = damping_dumin * 60
t = damping_t * 60

[slope, dslope, yitc, dyitc] = linreg(x=t, y=ln(f), dy=df/f)
D = -1.0 * slope
dD = dslope
f0 = exp(yitc)
df0 = exp(yitc) * dyitc

_f = f0 * exp(-D * t)
_df = (f0+df0) * exp(-(D+dD) * t)

dampingplot = plot(title='Calculation for damping constant', xlabel='t / s', ylabel='f / Hz', fig=1, scale='linlog')
dampingplot.plt.grid(False)
dampingplot.plotdata(t, f, df)
dampingplot.plotfunc(t, _f, label='line of best fit')
dampingplot.plotfunc(t, _df, label='line of uncertanty')

showfigs()
