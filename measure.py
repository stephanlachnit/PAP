### measure Python 3 libraby version 1.3.1
import math as m
import matplotlib.pyplot as plt
# todo: use numpy

# Settings
linreg_change = 0.00001 # min relative change per step to end linear regression

# Variables for export
sqrt = m.sqrt
exp = m.exp
ln = m.log
lg = m.log10
sin = m.sin
cos = m.cos
tan = m.tan
arctan = m.atan
pi = m.pi
euler_e = m.e
c = 2.99792458e8 # Speed of light
h = 6.626070040e-34 # Planck's constant
dh = 8.1e-42 # Uncertanty of Planck's constant
e = 1.6021766208e-19 # Elementary charge
de = 9.8e-28 # Uncertanty of the elementary charge
T0 = 273.15 # Zero Celsius in Kelvin

def mean_value(x):
  s = 0.0
  for i in range(len(x)):
    s += x[i]
  return s / len(x)

def std_dev_e(x):
  qs = 0.0
  for i in range(len(x)):
    qs += (x[i] - mean_value(x))**2
  return m.sqrt(qs / (len(x) - 1))

def std_dev_m(x):
  return std_dev_e(x) / m.sqrt(len(x))

def signval(val, err=0.0):
  if (err == 0.0):
    return "{:g}".format(val)
  if ("{:.1e}".format(err)[0] == "9" and "{:.0e}".format(err)[0] == "1"):
    err = float("{:.0e}".format(err))
  firstdigit = int("{:.1e}".format(err)[0])
  if (firstdigit <= 2):
    round2 = 1
    errstr = "{:.1e}".format(err)
  else:
    round2 = 0
    errstr = "{:.0e}".format(err)
  expdiff = int(m.floor(m.log10(abs(val))) - m.floor(m.log10(err)))
  if (expdiff < 0):
    sdigits = 0
    if (round2 != 1 or expdiff != -1):
      val = 0.0
  else:
    sdigits = expdiff + round2
  valstr = "{:.{digits}e}".format(val, digits=sdigits)
  return valstr + " ± " + errstr

def val(name, val, err=0.0):
  return name + ": " + signval(val, err)

def lst(name, val, err=[]):
  # todo: format data to make values align nicely, needs modifying of signval
  if (err == []):
    err = [0.0 for i in range(len(val))]
  tmp = name + ":"
  for i in range(len(val)):
    tmp +=  "\n " + signval(val[i], err[i])
  return tmp

def sig(name, val1, dVal1, val2, dVal2=0.0):
  nominator = abs(val1 - val2)
  denominator = m.sqrt(dVal1**2 + dVal2**2)
  if (nominator == 0.0):
    sigstr = "0"
  elif (denominator == 0.0):
    sigstr = "∞"
  else:
    sigma = nominator / denominator
    if (sigma < 0.95):
      digits = int(abs(m.floor(m.log10(sigma))))
    elif (sigma < 3.95):
      digits = 1
    else:
      digits = 0
    sigstr = "{:.{digits}f}".format(sigma, digits = digits)
  return name + ": " + sigstr + "σ"

def chi2(yo, dyo, ye, dye=[]):
  if (dye == []):
    dye = [0.0 for i in range(len(ye))]
  chi2 = 0.0
  for i in range(len(yo)):
    chi2 += (yo[i] - ye[i])**2 / (dyo[i]**2 + dye[i]**2)
  return chi2

def chi2_red(yo, dyo, ye, dye=[], dof=0):
  if (dof == 0):
    dof = len(ye)
  return chi2(yo, dyo, ye, dye) / dof

class plot:
  def __init__(self, title="", xlabel="", ylabel="", figure=1):
    self.figure = plt.figure(figure)
    self.legend = False
    plt.clf()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.ticklabel_format(style="sci", axis="both", scilimits=(-2,3))

  def plotdata(self, x, y, dy=[], dx=[], label=""):
    if (label != ""):
      self.legend = True
    if (dx == []):
      dx = [0.0 for i in range(len(x))]
    if (dy == []):
      dy = [0.0 for i in range(len(y))]
    plt.errorbar(x, y, dy, dx, label=label, fmt='o', markersize=3)
  
  def plotfunc(self, x, y, label=""):
    if (label != ""):
      self.legend = True
    plt.plot(x, y, label=label)

  def drawplot(self, show=True):
    if (show == True):
      if (self.legend == True):
        plt.legend()
    else:
      plt.close(self.figure)

  def saveplot(self, filename, fileformat="pdf"):
    plt.savefig(filename, format=fileformat)

  @staticmethod
  def showplots():
    plt.show()

def linreg(x, y, dy, dx=[], drawplot=False, title="", xlabel="", ylabel="", figure=1):
  def linreg_iter(x, y, dy):
    [s0, s1, s2, s3, s4] = [0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(len(x)):
      s0 += 1 / dy[i]**2
      s1 += x[i] / dy[i]**2
      s2 += y[i] / dy[i]**2
      s3 += x[i]**2 / dy[i]**2
      s4 += x[i] * y[i] / dy[i]**2
    eta = s0 * s3 - s1**2
    g = (s0 * s4 - s1 * s2) / eta
    dg = m.sqrt(s0 / eta)
    b = (s3 * s2 - s1 * s4) / eta
    db = m.sqrt(s3 / eta)
    return [g, dg, b, db]

  iter0 = linreg_iter(x, y, dy)
  result = []
  if (dx == []):
    result = iter0
  else:
    g = iter0[0]
    g_old = g * (1 - 2 * linreg_change)
    while (abs(1 - g_old / g) >= linreg_change):
      g_old = g
      dy = [m.sqrt((g * dx[i])**2 + dy[i]**2) for i in range(len(dy))]
      g = linreg_iter(x, y, dy)[0]
    result = linreg_iter(x, y, dy)
  if (drawplot == True):
    linregplot = plot(title=title, xlabel=xlabel, ylabel=ylabel, figure=figure)
    [g, dg, b, db] = result
    xn = len(x) - 1
    xint = [x[0] - dx[0], x[xn] + dx[xn]]
    yfit = [g * xint[i] + b for i in range(2)]
    yerr = [(g + dg) * xint[i] + (b - db) for i in range(2)]
    linregplot.plotdata(x, y, dy, dx)
    linregplot.plotfunc(xint, yfit, label="line of best fit")
    linregplot.plotfunc(xint, yerr, label="line of uncertainty")
    linregplot.drawplot()
    result.append(linregplot)
  return result

def lin_yerr(x, dx, y, dy):
  g = linreg(x, y, dx, dy)
  new_dy = [m.sqrt(dy[i]**2 + (g * dx[i])**2) for i in range(len(dy))]
  return new_dy
