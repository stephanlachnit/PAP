import math as m

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

def reg_grad(x, y, dy):
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

def reg_grad_err(x, y, dy):
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

def pv(name, val):
  print()
  print(name + " = " + str(val))

def pve(name, val, err):
  print()
  sigerr = float("{0:.1e}".format(err))
  sigval = round(val, 1 - int(m.floor(m.log10(sigerr))))
  print(name + " = " + str(sigval) + " ± " + str(sigerr))

def ple(name, val, err):
  print()
  print(name + ":")
  for i in range(len(val)):
    sigerr = float("{0:.1e}".format(err[i]))
    sigval = round(val[i], 1 - int(m.floor(m.log10(sigerr))))
    print(str(sigval) + " ± " + str(sigerr))
