# measure version 1.1

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
  x2 = 3.44
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

# Linear regression
blinreg = False
bplotlinreg = True
if (blinreg == True):
  from measure import linreg
  print()
  print("Linear regression")
  n = 10
  x = [i for i in range(n)]
  dx = [0.1 for i in range(n)]
  y = [2 + i**1.1 for i in range(n)]
  dy = [0.5 for i in range(n)]
  res = linreg(x, y, dy, dx, plot=bplotlinreg, title="measure test", xlabel="xlabel", ylabel="ylabel")
  [slope, slope_err, y_itc, y_itc_err] = res

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
