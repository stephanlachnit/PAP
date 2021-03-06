### measure libraby version 1.9.2s
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2 as _chi2
from numpy import loadtxt

# Settings
linreg_change = 0.00001 # min relative change per step to end linear regression
minfloat = 1e-80 # replaces zeros in linreg
linspace_res = 2000 # resolution for linspace

# Variables for export
sqrt = np.sqrt
exp = np.exp
ln = np.log
log10 = np.log10
sin = np.sin
cos = np.cos
tan = np.tan
arcsin = np.arcsin
arccos = np.arccos
arctan = np.arctan
pi = np.pi
euler_e = np.e
deg_to_rad = pi / 180
rad_to_deg = 180 / pi

c = 2.99792458e8 # Speed of light
h = 6.62607015e-34 # Planck's constant
e = 1.602176634e-19 # Elementary charge
kB = 1.380649e-23 # Boltzmann constant
NA = 6.02214076e23 # Avogadro constant
T0 = 273.15 # Zero Celsius in Kelvin
p0 = 101325 # NIST standard pressure
g = 9.80984 # Gravitanional acceleration in Heidelberg 
dg = 2e-5 # Uncertainty of the gravitational acceleration

def npfarray(x):
  return np.array(x, dtype='float')

def nplinspace(start, stop):
  return np.linspace(start,stop,num=linspace_res,endpoint=True)

def spcurvefit(fitfunc, xdata, ydata, yerr=None, p0=None):
  p_opt,p_cov = curve_fit(fitfunc, xdata, ydata, p0=p0, sigma=yerr)
  p_err = sqrt(np.diag(p_cov))
  return (p_opt, p_err)

def mv(val):
  """
  Parameters

  val: npfarray (or list);
  ----------
  Returns

  float; mean value of val;
  """
  if (type(val) == list):
    val = npfarray(val)
  s = np.sum(val)
  return s / len(val)

def dsto(val, ddof=1):
  """
  Parameters

  val: npfarray (or list);
  ddof: float;
  ----------
  Returns

  float; standard devation with dof = len(val) - ddof;
  """
  if (type(val) == list):
    val = npfarray(val)
  _mv = mv(val)
  s = np.sum((val - _mv)**2)
  return sqrt(s / (len(val) - ddof))

def dsto_mv(val, ddof=1):
  """
  Parameters

  val: npfarray (or list);
  ddof: float;
  ----------
  Returns

  float; standard devation of the mean value with dof = len(val) - ddof;
  """
  return dsto(val,ddof=ddof) / sqrt(len(val))

def dsys_mv(err):
  """
  Parameters

  err: npfarray (or list); systematic errors;
  ----------
  Returns

  float; systematic error of the mean value;
  """
  if (type(err) == list):
    err = npfarray(err)
  return sqrt(np.sum(err**2)) / len(err)

def dtot(dsys, dsto):
  """
  Parameters

  dsys: float; systematic error;
  dsto: float; stochastic error;
  ----------
  Returns

  float; total error;
  """
  return sqrt(dsys**2 + dsto**2)

def dtot_mv(val, err):
  """
  Parameters

  val: npfarray (or list);
  err: npfarray (or list); systematic errors;
  ----------
  Returns

  float; total error of the mean value;
  """
  return dtot(dsys_mv(err),dsto_mv(val))

def signval(val, err=0.0):
  if (err == 0.0):
    return ['{:g}'.format(val), '']
  errstr = '{:.1e}'.format(err)
  err = float(errstr)
  expdiff = int(np.floor(np.log10(abs(val))) - np.floor(np.log10(err)))
  sdigits = expdiff + 1
  if (expdiff < -2):
    val = 0.0
    sdigits = 0
  elif (expdiff == -2):
    if (int('{:e}'.format(abs(val))[0]) >= 5):
      sig = val / abs(val)
      val = sig * 10**int(1+np.floor(np.log10(abs(val))))
    else:
      val = 0.0
    sdigits = 0
  elif (expdiff == -1):
    if (int(np.floor(np.log10(abs(float('{:.0e}'.format(val))))) - np.floor(np.log10(err))) == -1):
      val = val
      sdigits = 0
    else:
      val = float('{:.0e}'.format(val))
      sdigits = 1
  else:
    if (int(np.floor(np.log10(abs(float('{:.{digits}e}'.format(val, digits=sdigits))))) - np.floor(np.log10(err))) != expdiff):
      val = float('{:.{digits}e}'.format(val, digits=sdigits))
      sdigits += 1
  valstr = '{:.{digits}e}'.format(val, digits=sdigits)
  return [valstr, errstr]

def val(val, err=0.0, name='', percerr=False):
  """
  Parameters

  val: float;
  err: float; uncertainty of val;
  name: string; name of val;
  percerr: bool; print error in percent
  ----------
  Returns

  string; format: "name = val ± err" with two significant digits;
  """
  out = ''
  if (name != ''):
    out += name + ' = '
  tmp = signval(val, err)
  out += tmp[0]
  if (tmp[1] != ''):
    out += ' ± ' + tmp[1]
    if percerr:
      out += ' ({:.2g}%)'.format(100 * err / val)
  return out

def lst(val, err=[], name=''):
  """
  Parameters

  val: npfarray;
  err: npfarray; uncertainties of val;
  name: string; name of the list;
  ----------
  Returns

  array of strings; format "val[i] ± err[i]" with two significant digits, name is the first item;
  """
  if (err == []):
    err = [0.0 for i in range(len(val))]
  N = len(val)
  valmaxlen = 0
  errmaxlen = 0
  hasneg = False
  for i in range(N):
    tmp = signval(val[i], err[i])
    if (len(tmp[0]) > valmaxlen):
      valmaxlen = len(tmp[0])
    if (len(tmp[1]) > errmaxlen):
      errmaxlen = len(tmp[1])
    if (float(tmp[0]) < 0):
      hasneg = True
  out = []
  if (name != ''):
    erraddlen = 0
    if (errmaxlen != 0):
      erraddlen = 3
    adjust = int(np.floor((valmaxlen + errmaxlen + erraddlen + len(name))/2))
    out.append(name.rjust(adjust))
  for i in range(len(val)):
    tmp = signval(val[i], err[i])
    if hasneg:
      if (float(tmp[0]) >= 0):
        tmp[0] = '+' + tmp[0]
    tmp2 = tmp[0].ljust(valmaxlen)
    if (tmp[1] != ''):
      tmp2 += ' ± ' + tmp[1].ljust(errmaxlen)
    elif (errmaxlen != 0):
      tmp2 += ''.ljust(errmaxlen + 3)
    out.append(tmp2)
  return out

def tbl(lists):
  """
  Parameters

  lists: array of arrays of strings; each list (stringarray) will turn into a column;
  ----------
  Returns

  string; format: table;
  """
  M = len(lists[0])
  N = len(lists)
  out = ''
  lens = [int(npfarray([len(lists[i][j]) for j in range(M)]).max()) for i in range(N)]
  for j in range(M):
    for i in range(N):
      suffix = ' | '
      if (i == N-1):
        suffix = '\n'
        if (j == M-1):
          suffix = ''
      out += lists[i][j].ljust(lens[i]) + suffix
  return out

def dev(val1, err1, val2, err2=0.0, name='', perc=False):
  def get_sig(nominator,denominator):
    if (nominator == 0.0):
      sigstr = '0'
    elif (denominator == 0.0):
      sigstr = '∞ '
    else:
      sigma = nominator / denominator
      if (sigma < 0.95):
        digits = int(abs(np.floor(np.log10(sigma))))
      elif (sigma < 3.95):
        digits = 1
      else:
        digits = 0
      sigstr = '{:.{digits}f}'.format(sigma,digits=digits)
    sigstr += 'σ'
    return sigstr
  def get_perc(val1,val2,pformat='{:.1f}'):
    percval = abs(val1 - val2) / val2 * 100
    percstr = pformat.format(percval) + '%'
    return percstr
  out = None
  array = False
  if type(val1) is list:
    val1 = npfarray(val1)
    array = True
  elif type(val1) is np.ndarray:
    array = True
  if type(val2) is list:
    val2 = npfarray(val2)
    array = True
  elif type(val2) is np.ndarray:
    array = True
  if type(err1) is list:
    err1 = npfarray(err1)
    array = True
  elif type(err1) is np.ndarray:
    array = True
  if type(err2) is list:
    err2 = npfarray(err2)
    array = True
  elif type(err2) is np.ndarray:
    array = True
  nominator = abs(val1 - val2)
  denominator = np.sqrt(err1**2 + err2**2)
  if array:
    out = []
    N = len(val1)
    if type(denominator) is not np.ndarray:
      denominator = npfarray([denominator for i in range(N)])
    if perc:
      if type(val2) is not np.ndarray:
        val2 = npfarray([val2 for i in range(N)])
      tmp = []
      tmp2 = []
      sigmaxlen = 0
      percmaxlen = 0
      for i in range(N):
        tmp.append(get_sig(nominator[i],denominator[i]))
        siglen = len(tmp[i])
        if (siglen > sigmaxlen):
          sigmaxlen = siglen
        tmp2.append(get_perc(val1[i],val2[i]))
        perclen = len(tmp2[i])
        if (perclen > percmaxlen):
          percmaxlen = perclen
      if (name != ''):
        adjust = int(np.floor((sigmaxlen/2 + 2 + len(name))))
        out.append(name.rjust(adjust))
      for i in range(N):
        out.append(tmp[i].rjust(sigmaxlen) + ' | ' + tmp2[i].rjust(percmaxlen))
    else:
      if (name != ''):
        out.append(name)
      for i in range(N):
        out.append(get_sig(nominator[i],denominator[i]))
  else:
    out = ''
    prefix = ''
    if (name != ''):
      prefix = name + ': '
    out += prefix + get_sig(nominator,denominator)
    if perc:
      out += ' ; ' + get_perc(val1,val2,pformat='{:.2g}')
  return out

class pltext:
  @staticmethod
  def initplot(num=0, title='', xlabel='', ylabel='', scale='linlin'):
    fig = plt.figure(num)
    plt.title(title, fontsize='14')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, which='both')
    fig.set_size_inches(11.69,8.27)
    if (scale == 'linlin'):
      plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,3))
    elif (scale == 'linlog'):
      plt.yscale('log')
      plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,3))
    elif (scale == 'loglin'):
      plt.xscale('log')
      plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,3))
    elif (scale == 'loglog'):
      plt.yscale('log')
      plt.xscale('log')

  @staticmethod
  def plotdata(x, y, dy=None, dx=None, label='', caps=True, connect=False, color=None):
    capsize = 0
    if caps:
      capsize = 3
    plot = plt.errorbar(x=x,y=y,yerr=dy,xerr=dx,label=label,color=color,fmt='o',markersize=3,capsize=capsize)
    if (color == None):
      color = plot[0].get_color()
    if (connect == True):
      plt.plot(x, y, color=color)
    return color
  
  @staticmethod
  def set_layout(xlim=None, ylim=None, legend=True):
    if (xlim != None):
      plt.xlim(xlim)
    if (ylim != None):
      plt.ylim(ylim)
    if legend:
      plt.legend()
    plt.tight_layout()

class chi2stat:
  @staticmethod
  def chi2(val1, err1, val2, err2=0.0):
    return np.sum((val1 - val2)**2 / (err1**2 + err2**2))

  @staticmethod
  def chi2_red(chi2, len, ddof=0):
    return chi2 / (len - ddof)

  @staticmethod
  def fit_prob(chi2, len, ddof=0):
    return 1 - _chi2.cdf(chi2, len - ddof)

def linreg(x, y, dy, dx=None, plot=False, prange=None, uncertdraw='slope_dec', label=''):
  def linreg_iter(x, y, dy):
    [s0, s1, s2, s3, s4] = [0.0, 0.0, 0.0, 0.0, 0.0]
    for i in range(len(x)):
      if (dy[i] == 0.0):
        dy[i] = minfloat
      s0 += dy[i]**-2
      s1 += x[i] * dy[i]**-2
      s2 += y[i] * dy[i]**-2
      s3 += x[i]**2 * dy[i]**-2
      s4 += x[i] * y[i] * dy[i]**-2
    eta = s0 * s3 - s1**2
    g = (s0 * s4 - s1 * s2) / eta
    dg = np.sqrt(s0 / eta)
    b = (s3 * s2 - s1 * s4) / eta
    db = np.sqrt(s3 / eta)
    return [g, dg, b, db]
  iter0 = linreg_iter(x, y, dy)
  result = []
  _dx = dx
  if (dx == None):
    _dx = np.zeros(len(x))
    result = iter0
  else:
    g = iter0[0]
    g_old = g * (1 - 2 * linreg_change)
    _dy = dy
    while (abs(1 - g_old / g) >= linreg_change):
      g_old = g
      _dy = np.sqrt((g * dx)**2 + _dy**2)
      g = linreg_iter(x, y, _dy)[0]
    result = linreg_iter(x, y, _dy)
  if plot:
    if (prange == None):
      xmin_index = np.argmin(x - _dx)
      xmax_index = np.argmax(x + _dx)
      xmin = x[xmin_index] - _dx[xmin_index]
      xmax = x[xmax_index] + _dx[xmax_index]
      xint = nplinspace(xmin, xmax)
    else:
      xint = nplinspace(*prange)
    prefix = ''
    color = None
    if (label != ''):
      prefix = label + ': '
      color = pltext.plotdata(x=x, y=y, dy=dy, dx=dx, label=prefix+'Measurements')
    [g, dg, b, db] = result
    yfit = g * xint + b
    color = plt.plot(xint, yfit, marker='', label=prefix+'Fit', color=color)[0].get_color()
    if (uncertdraw == 'slope_inc'):
      yerr = (g + dg) * xint + (b - db)
      plt.plot(xint, yerr, marker='', linestyle='dashed', label=prefix+'Uncertainty', color=color)
    elif (uncertdraw == 'slope_dec'):
      yerr = (g - dg) * xint + (b + db)
      plt.plot(xint, yerr, marker='', linestyle='dashed', label=prefix+'Uncertainty', color=color)
  return result

def expreg(x, y, dy, dx=None, plot=False, prange=None, uncertdraw='expon_inc', label=''):
  expo,dexpo,_yitc,_dyitc = linreg(x,np.log(y),dy/y,dx,plot=False)
  yitc = exp(_yitc)
  dyitc = yitc * _dyitc
  result = [expo,dexpo,yitc,dyitc]
  if (plot):
    _dx = dx
    if (dx == None):
      _dx = np.zeros(len(x))
    if (prange == None):
      xmin_index = np.argmin(x - _dx)
      xmax_index = np.argmax(x + _dx)
      xmin = x[xmin_index] - _dx[xmin_index]
      xmax = x[xmax_index] + _dx[xmax_index]
      xint = nplinspace(xmin, xmax)
    else:
      xint = nplinspace(*prange)
    prefix = ''
    color = None
    if (label != ''):
      prefix = label + ': '
      color = pltext.plotdata(x=x, y=y, dy=dy, dx=dx, label=prefix+'Measurements')
    yfit = yitc * exp(expo*xint)
    yerr = (yitc-dyitc) * exp((expo+dexpo) * xint)
    color = plt.plot(xint, yfit, marker='', label=prefix+'Fit', color=color)[0].get_color()
    if (uncertdraw == 'expon_inc'):
      yerr = (yitc-dyitc) * exp((expo+dexpo) * xint)
      plt.plot(xint, yerr, marker='', linestyle='dashed', label=prefix+'Uncertainty', color=color)
    elif (uncertdraw == 'expon_dec'):
      yerr = (yitc+dyitc) * exp((expo-dexpo) * xint)
      plt.plot(xint, yerr, marker='', linestyle='dashed', label=prefix+'Uncertainty', color=color)    
  return result
