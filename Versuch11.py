import math

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
  return math.sqrt(s0 / eta)

"""
m = [50, 100, 150, 200, 250]
t = [0.92, 1.25, 1.47, 1.66, 1.85]

x = m
y = [x * x for x in t]
dy = [0.028, 0.007, 0.04, 0.023, 0.04]

g = reg_grad(x, y, dy)
dg = reg_grad_err(x, y, dy)

print("g = " + str(g))
print("dg = " + str(dg))
print("d = " + str(4 * math.pi * math.pi / g))
print("dd = " + str(4 * math.pi * math.pi / (g * g) * dg))
"""

d = 3192.128167101808
dd = 43.507500508062634

m = [0, 50, 100, 150, 200, 250]
s = [45, 207, 372, 537, 699, 858]

x = m
y = s
dy = [1, 1, 1, 1, 1, 1]

a = reg_grad(x, y, dy)
da = reg_grad_err(x, y, dy)

print("a = " + str(a))
print("da = " + str(da))
print("g = " + str(d * a))
print("dg = " + str(math.sqrt(a * a * dd * dd + + d * d * da * da)))