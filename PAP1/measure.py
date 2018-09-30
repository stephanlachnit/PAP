import math as m
import matplotlib.pyplot as plt

p = 0.00001 # percent change to stop iterations in linear regression
n = 1000 # number of iteration for plotting functions

def sqrt(x):
  return m.sqrt(x)

def pi():
  return m.pi

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

def reg_itc_y(x, y, dy):
  s0 = 0
  s1 = 0
  s2 = 0
  s3 = 0
  s4 = 0
  for i in range(len(x)):
    s0 += 1 / (dy[i] * dy[i])
    s1 += x[i] / (dy[i] * dy[i])
    s2 += y[i] / (dy[i] * dy[i])
    s3 += x[i] * x[i] / (dy[i] * dy[i])
    s4 += x[i] * y[i] / (dy[i] * dy[i])
  eta = s0 * s3 - s1 * s1
  return (s3 * s2 - s1 * s4) / eta

def reg_itc_err_y(x, y, dy):
  s0 = 0
  s1 = 0
  s2 = 0
  s3 = 0
  s4 = 0
  for i in range(len(x)):
    s0 += 1 / (dy[i] * dy[i])
    s1 += x[i] / (dy[i] * dy[i])
    s2 += y[i] / (dy[i] * dy[i])
    s3 += x[i] * x[i] / (dy[i] * dy[i])
    s4 += x[i] * y[i] / (dy[i] * dy[i])
  eta = s0 * s3 - s1 * s1
  return m.sqrt(s3 / eta)

def reg_grad_y(x, y, dy):
  s0 = 0
  s1 = 0
  s2 = 0
  s3 = 0
  s4 = 0
  for i in range(len(x)):
    s0 += 1 / (dy[i] * dy[i])
    s1 += x[i] / (dy[i] * dy[i])
    s2 += y[i] / (dy[i] * dy[i])
    s3 += x[i] * x[i] / (dy[i] * dy[i])
    s4 += x[i] * y[i] / (dy[i] * dy[i])
  eta = s0 * s3 - s1 * s1
  return (s0 * s4 - s1 * s2) / eta

def reg_grad_err_y(x, y, dy):
  s0 = 0
  s1 = 0
  s2 = 0
  s3 = 0
  s4 = 0
  for i in range(len(x)):
    s0 += 1 / (dy[i] * dy[i])
    s1 += x[i] / (dy[i] * dy[i])
    s2 += y[i] / (dy[i] * dy[i])
    s3 += x[i] * x[i] / (dy[i] * dy[i])
    s4 += x[i] * y[i] / (dy[i] * dy[i])
  eta = s0 * s3 - s1 * s1
  return m.sqrt(s0 / eta)

def reg_itc(x, y, dx, dy):
  g = reg_grad(x, y, dx, dy, False)
  newDy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
  return reg_itc_y(x, y, newDy)

def reg_itc_err(x, y, dx, dy):
  g = reg_grad(x, y, dx, dy, False)
  newDy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
  return reg_itc_err_y(x, y, newDy)

def reg_grad(x, y, dx, dy, print_g=False):
  g = reg_grad_y(x, y, dy)
  if print_g:
    print("g = " + str(g))
  gOld = 2 * g
  while abs(1 - gOld / g) > p:
    gOld = g
    newDy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
    g = reg_grad_y(x, y, newDy)
    if print_g:
      print("g = " + str(g))
  return g

def reg_grad_err(x, y, dx, dy):
  g = reg_grad(x, y, dx, dy, False)
  newDy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
  return reg_grad_err_y(x, y, newDy)

def chi2(ve, vo, dv):
  ch2 = 0
  for i in range(len(x)):
    ch2 += (vo - ve)**2 / dy[i]**2
  return ch2 / len(x);

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

def plot(title, xlabel, ylabel, xval, xerr, yval, yerr):
  xBeg = xval[0]
  xEnd = xval[len(xval)-1]

  g = reg_grad(xval, yval, xerr, yerr)
  dg = reg_grad_err(xval, yval, xerr, yerr)
  b = reg_itc(xval, yval, xerr, yerr)
  db = reg_itc_err(xval, yval, xerr, yerr)

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
  plt.plot(x, y0, label="line of fit")
  plt.plot(x, y1, label="line of error")
  plt.errorbar(xval, yval, yerr, xerr, fmt='x', capsize=2.0)
  plt.legend()
  plt.show()
  return [g, dg, b, db]
