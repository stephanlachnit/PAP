# measure version 1.6

# Significant values
bval = False
if (bval == True):
  print()
  print('Significant values')
  from measure import val
  x0 = 3.357789123e4
  print(val('no error', x0))
  x1 = 760.0e-9
  d1 = 0.11e-9
  print(val('first error digit <= 2', x1, d1))
  x2 = 531.2e3
  d2 = 3.1e3
  print(val('first error digit >= 3', x2, d2))

# Tables
btable = False
if (btable == True):
  print()
  print('Tables')
  from measure import tbl
  data = [ [1,2,3], [2,4,6] ]
  names = ['a', 'b']
  print(tbl(names, data))

# Statistical evaluation
bstat = False
if (bstat == True):
  print()
  print('Statistical evaluation')
  from measure import val,mean_value,std_dev_e,std_dev_m,chi2, chi2_red
  x = [2.33, 2.37, 2.21, 2.32, 2.17, 2.44, 2.11]
  x0 = x[2]
  xm = mean_value(x)
  dx = std_dev_e(x)
  dxm = std_dev_m(x)
  print(val('Mean value', xm, dxm))
  print(val('Some value', x0, dx))
  xobserved = [1.1, 1.9, 3.2, 3.8, 5.2]
  dxo = [0.1, 0.3, 0.2, 0.2, 0.3]
  xexpected = [1.0, 2.0, 3.0, 4.0, 5.0]
  dxt = [1e-2, 2e-2, 3e-2, 4e-2, 5e-2]
  _chi2 = chi2(xobserved, dxo, xexpected, dye=dxt)
  number_of_measurements = len(xobserved)
  dependency_factor = 2
  degrees_of_freedom = number_of_measurements - dependency_factor
  _chi2_red = chi2_red(xobserved, dxo, xexpected, dye=dxt, dof=degrees_of_freedom)
  print(val('𝜒²', _chi2))
  print(val('𝜒² reduced', _chi2_red))

# Deviation between values
bsig = False
if (bsig == True):
  print()
  print('Deviation between values')
  from measure import sig
  x1 = 3.4
  d1 = 0.1
  x2 = 3.404
  print(sig('small deviation', x1, d1, x2))
  x1 = 9.7e3
  d1 = 1.7e3
  x2 = 4.2e3
  d2 = 1.2e-1
  print(sig('big deviation', x1, d1, x2, dVal2=d2))

# Constants
bconstants = False
if (bconstants == True):
  print()
  print('Constants')
  from measure import val,pi,euler_e as e
  print(val('Pi', pi))
  print(val('e', e))
  from measure import c,e,de,h,dh,T0,g,dg
  print(val('Speed of light', c))
  print(val('Elementary charge', e, de))
  print(val('Planck\'s constant', h, dh))
  print(val('Zero Celsiuis in Kelvin', T0))
  print(val('Gravitational acceleration in Heidelberg', g, dg))

# Plot data
bplotdata = False
if (bplotdata == True):
  print()
  print('Plot data')
  from measure import plot
  x = [1, 2, 3]
  dx = [0.1, 0.2, 0.3]
  y = [8, 2, 1]
  dy = [2, 1, 0.5]
  dataplot = plot(title='data plot', xlabel='x', ylabel='y', fig=1)
  dataplot.plotdata(x, y, dy, dx, label='samples')
  plot.showfigs()

# Plot functions
bplotfunc = False
if (bplotfunc == True):
  print()
  print('Plot functions')
  from measure import plot
  x = [0.1 * i for i in range(100)]
  y = [i**2 for i in range(100)]
  funcplot = plot(title='function plot', xlabel='x', ylabel='y', fig=2)
  funcplot.plotfunc(x, y, label='f(x)=x^2')
  plot.showfigs()

# Linear regression w/o plot
blinreg = False
if (blinreg == True):
  print()
  print('Linear regression w/o plot')
  from measure import linreg,val
  n = 10
  x = [i for i in range(n)]
  dx = [0.1 for i in range(n)]
  y = [2 + i**1.1 for i in range(n)]
  dy = [0.5 for i in range(n)]
  [slope, slope_err, y_itc, y_itc_err] = linreg(x, y, dy, dx)
  print(val('slope', slope, slope_err))
  print(val('y_intercept', y_itc, y_itc_err))

# Linear regression w/ plot
blinregplot = False
if (blinregplot == True):
  print()
  print('Linear regression w/ plot')
  from measure import linreg,plot
  n = 10
  x = [i for i in range(n)]
  dx = [0.1 for i in range(n)]
  y = [2 + i**0.9 for i in range(n)]
  dy = [0.5 for i in range(n)]
  lrplot = plot(title='Linear regression', xlabel='x', ylabel='y', fig=3)
  [g, dg, b, db, linregplot] = linreg(x, y, dy, dx, lrplot=lrplot)
  plot.showfigs()
