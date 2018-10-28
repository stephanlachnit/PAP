# measure version 1.3

# Significant values
bval = True
if (bval == True):
  from measure import val
  print()
  print("Significant values")
  x0 = 3.357789123e4
  print(val("no error", x0))
  x1 = 760.0e-9
  d1 = 0.11e-9
  print(val("first error digit <= 2", x1, d1))
  x2 = 531.2e3
  d2 = 3.1e3
  print(val("first error digit >= 3", x2, d2))

# Statistical evaluation
bstat = False
if (bstat == True):
  from measure import val
  print()
  print("Statistical evaluation")
  from measure import mean_value,std_dev_e,std_dev_m
  x = [2.33, 2.37, 2.21, 2.32, 2.17, 2.44, 2.11]
  x0 = x[2]
  xm = mean_value(x)
  dx = std_dev_e(x)
  dxm = std_dev_m(x)
  print(val("Mean value", xm, dxm))
  print(val("Some value", x0, dx))
  from measure import chi2, chi2_red
  xobserved = [1.1, 1.9, 3.2, 3.8, 5.2]
  dxo = [0.1, 0.3, 0.2, 0.2, 0.3]
  xexpected = [1.0, 2.0, 3.0, 4.0, 5.0]
  dxt = [1e-2, 2e-2, 3e-2, 4e-2, 5e-2]
  _chi2 = chi2(xobserved, dxo, xexpected, dye=dxt)
  number_of_measurements = len(xobserved)
  dependency_factor = 2
  degrees_of_freedom = number_of_measurements - dependency_factor
  _chi2_red = chi2_red(xobserved, dxo, xexpected, dye=dxt, dof=degrees_of_freedom)
  print(val("𝜒²", _chi2))
  print(val("𝜒² reduced", _chi2_red))

# Deviation between values
bsig = False
if (bsig == True):
  from measure import sig
  print()
  print("Deviation between values")
  x1 = 3.4
  d1 = 0.1
  x2 = 3.404
  print(sig("small deviation", x1, d1, x2))
  x1 = 9.7e3
  d1 = 1.7e3
  x2 = 4.2e3
  d2 = 1.2e-1
  print(sig("big deviation", x1, d1, x2, dVal2=d2))

# Constants
bconstants = False
if (bconstants == True):
  from measure import val
  print()
  print("Constants")
  from measure import pi
  print(val("Pi", pi))
  from measure import euler_e as e
  print(val("e", e))
  from measure import c
  print(val("Speed of light", c))
  from measure import e,de
  print(val("Elementary charge", e, de))
  from measure import h,dh
  print(val("Plancl's constant", h, dh))

# Plot data
bplotdata = False
if (bplotdata == True):
  from measure import plot
  print()
  print("Plot data")
  x = [1, 2, 3]
  dx = [0.1, 0.2, 0.3]
  y = [8, 2, 1]
  dy = [2, 1, 0.5]
  dataplot = plot(title="data plot", xlabel="x", ylabel="y", figure=1)
  dataplot.plotdata(x, y, dy, dx, label="samples")
  dataplot.drawplot()
  plot.showplots()

# Plot functions
bplotfunc = False
if (bplotdata == True):
  from measure import plot
  print()
  print("Plot functions")
  x = [0.1 * i for i in range(100)]
  y = [i**2 for i in range(100)]
  funcplot = plot(title="function plot", xlabel="x", ylabel="y", figure=1)
  funcplot.plotfunc(x, y, label="f(x)=x^2")
  funcplot.drawplot()
  plot.showplots()

# Linear regression w/o plot
blinreg = False
if (blinreg == True):
  from measure import linreg,val
  print()
  print("Linear regression w/o plot")
  n = 10
  x = [i for i in range(n)]
  dx = [0.1 for i in range(n)]
  y = [2 + i**1.1 for i in range(n)]
  dy = [0.5 for i in range(n)]
  [slope, slope_err, y_itc, y_itc_err] = linreg(x, y, dy, dx)
  print(val("slope", slope, slope_err))
  print(val("y_intercept", y_itc, y_itc_err))

# Linear regression w/ plot
blinregplot = False
if (blinregplot == True):
  from measure import linreg,plot
  print()
  print("Linear regression w/ plot")
  n = 10
  x = [i for i in range(n)]
  dx = [0.1 for i in range(n)]
  y = [2 + i**0.9 for i in range(n)]
  dy = [0.5 for i in range(n)]
  [g, dg, b, db, linregplot] = linreg(x, y, dy, dx, drawplot=True, title="Linear regression")
  plot.showplots()
