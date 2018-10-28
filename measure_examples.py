# measure version 1.3

# Significant values
bval = False
if (bval == True):
  from measure import val
  print()
  print("Significant values")
  x1 = 760.0e-9
  d1 = 0.11e-9
  print(val("first error digit <= 2", x1, d1))
  x2 = 531.2e3
  d2 = 3.1e3
  print(val("first error digit >= 3", x2, d2))

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
  print(sig("big deviation", x1, d1, x2, d2))

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

# Plot functions
bplotfunc = False
if (bplotdata == True):
  from measure import plot
  print()
  print("Plot functions")

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
