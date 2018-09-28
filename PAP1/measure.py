import math as m

p = 0.001

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
  qs = 0
  for i in range(len(x)):
    qs += (x[i] - mean_value(x))**2
  return m.sqrt(qs / ((len(x) - 1) * len(x)))

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

def reg_grad(x, y, dx, dy):
  g = reg_grad_y(x, y, dy)
  print("g = " + str(g))
  gOld = 2 * g
  while abs(1 - gOld / g) > p:
    gOld = g  
    g = reg_grad_y(x, y, m.sqrt((g * dx)**2 + dy**2))
    print("g = " + str(g))
  return g

def reg_grad_err(x, y, dx, dy):
  g = reg_grad(x, y, dx, dy)
  return reg_grad_err_y(x, y, m.sqrt((g * dx)**2 + dy**2))

def sigerr(err):
  if ("{0:.1e}".format(err)[0] == "1" or "{0:.1e}".format(err)[0] == "2"):
    return float("{0:.1e}".format(err))
  else:
    return float("{0:.0e}".format(err))

def sigval(val, err):
  return round(val, 1 - int(m.floor(m.log10(sigerr(err)))))

def pv(name, val):
  print()
  print(name + " = " + str(val))

def pve(name, val, err):
  print()
  print(name + " = " + str(sigval(val, err)) + " ± " + str(sigerr(err)))

def pl(name, val):
  print()
  print(name + ":")
  for i in range(len(val)):
    print(val[i])

def ple(name, val, err):
  print()
  print(name + ":")
  for i in range(len(val)):
    print(str(sigval(val[i], err[i])) + " ± " + str(sigerr(err[i])))
