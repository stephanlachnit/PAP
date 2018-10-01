import math as m
import matplotlib.pyplot as plt

p = 0.00001 # percent change to stop iterations in linear regression
n = 1000 # number of iteration for plotting functions

pi = m.pi
sqrt = m.sqrt
ln = m.log

def mean_value(x):
  s = 0
  for i in range(len(x)):
    s += x[i]
  return s / len(x)

def std_dev_e(x):
  qs = 0
  for i in range(len(x)):
    qs += (x[i] - mean_value(x))**2
  return m.sqrt(qs / (len(x) - 1))

def std_dev_m(x):
  return std_dev_e(x) / sqrt(len(x))

def linreg_iter(x, y, dy):
  s0 = 0
  s1 = 0
  s2 = 0
  s3 = 0
  s4 = 0
  for i in range(len(x)):
    s0 += 1 / dy[i]**2
    s1 += x[i] / dy[i]**2
    s2 += y[i] / dy[i]**2
    s3 += x[i]**2 / dy[i]**2
    s4 += x[i] * y[i] / dy[i]**2
  eta = s0 * s3 - s1**2
  slope = (s0 * s4 - s1 * s2) / eta
  slope_err = m.sqrt(s0 / eta)
  y_itc = (s3 * s2 - s1 * s4) / eta
  y_itc_err = m.sqrt(s3 / eta)
  return [slope, slope_err, y_itc, y_itc_err]

def linreg(x, y, dy, dx = []):
  iter0 = linreg_iter(x, y, dy)
  if (dx == []):
    return iter0
  else:
    g = iter0[0]
    g_old = g * (1 - 2 * p)
    while (abs(1 - g_old / g) >= p):
      g_old = g
      dy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
      g = linreg_iter(x, y, dy)[0]
    return linreg_iter(x, y, dy)

def lin_yerr(x, dx, y, dy):
  g = linreg(x, y, dx, dy)
  new_dy = [m.sqrt(dy[i]**2 + (g * dx[i])**2) for i in range(len(dy))]
  return new_dy

def chi2(yo, dyo, ye, dye = []):
  if (dye == []):
    dye = [0.0 for i in range(len(ye))]
  x2 = 0.0
  for i in range(len(yo)):
    x2 += (yo[i] - ye[i])**2 / (dyo[i]**2 + dye[i]**2)
  return x2 / len(ye)

def sigerr(err):
  if ("{0:.1e}".format(err)[0] == "1" or "{0:.1e}".format(err)[0] == "2"):
    return float("{0:.1e}".format(err))
  else:
    return float("{0:.0e}".format(err))

def sigval(val, err):
  _sigerr = sigerr(err)
  if (_sigerr == 0):
    return val
  else:
    if (int("{0:.0e}".format(err)[0]) < 3):
      round2 = 1
    else:
      round2 = 0
    return round(val, round2 - int(m.floor(m.log10(_sigerr))))

def pv(name, val, space=True):
  print(name + ": " + str(val))
  if space:
    print()

def pve(name, val, err, space=True):
  print(name + ": " + str(sigval(val, err)) + " ± " + str(sigerr(err)))
  if space:
    print()

def pl(name, val, space=True):
  print(name + ":")
  for i in range(len(val)):
    print(val[i])
  if space:
    print()

def ple(name, val, err, space=True):
  print(name + ":")
  for i in range(len(val)):
    print(str(sigval(val[i], err[i])) + " ± " + str(sigerr(err[i])))
  if space:
    print()

def ps(name, val1, val2, dVal1, dVal2, space=True):
  dVal = max([dVal1, dVal2])
  sigma = abs(val1 - val2) / dVal
  if sigma == 0:
    print(name + ": " + str(0.0) + "σ")
  else:
    print(name + ": " + str(sigval(sigma, pow(10, int(m.log10(sigma)) - 1))) + "σ")
  if space:
    print()

def plot(title, xlabel, ylabel, xval, yval, yerr, xerr = []):
  if (xerr == []):
    xerr = [0.0 for i in range(len(xval))]
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.plot(xval, yval)
  plt.errorbar(xval, yval, yerr, xerr, fmt='x', capsize=2.0)
  plt.legend()
  plt.show()

def plot_linreg(title, xlabel, ylabel, xval, yval, yerr, xerr = []):
  if (xerr == []):
    xerr = [0.0 for i in range(len(xval))]
  xBeg = xval[0]
  xEnd = xval[len(xval)-1]
  [g, dg, b, db] = linreg(xval, yval, yerr, xerr)
  x = []
  y0 = []
  y1 = []
  for i in range(n):
    x.append(xBeg + (xEnd - xBeg) * i / (n - 1))
    y0.append( b + g * x[i])
    y1.append((b - db) + (g + dg) * x[i])
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.plot(x, y0, label="line of best fit")
  plt.plot(x, y1, label="line of uncertanty")
  plt.errorbar(xval, yval, yerr, xerr, fmt='x', capsize=2.0)
  plt.legend()
  plt.show()
  return [g, dg, b, db]
