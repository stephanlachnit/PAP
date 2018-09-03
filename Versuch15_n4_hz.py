import measure.py

g = 

alpha = 
dalpha = 
r1 = 
dr1 = 
r2 = 
dr2 = 

t_m = [[], [], [], []]
s = [ , , , ]
ds = [ , , , ]

t = [mean_value(x) for x in t_m]
dt = [std_dev_m(x) for x in t_m]

t2 = [x**2 for x in t]
dt2 = [abs(2*t[i] * dt[i]) for i in range(len(dt))]

a_m = 2.0 * reg_grad(s, t2, dt2)
da_m = 2.0 * reg_grad_err(s, t2, dt2)


a = 2.0 * g * sin(alpha) / (3.0 - (r1 / r2)**2)
da = sqrt((2.0 * g * cos(alpha) / (3.0 - (r1 / r2)**2) * dalpha)**2
           + (4.0 * g * sin(alpha) / (3.0 - (r1 / r2)**2)**2 * r1 / (r2**2) * dr1)**2
           + (4.0 * g * sin(alpha) / (3.0 - (r1 / r2)**2)**2 * r1**2 / (r2**3) * dr2)**2)

print("Gemessene Beschleunigung: " + a_m)
print("Fehler der gemessenen Beschleunigung: " + da_m)
print("Aus Messdaten errechnete Beschleunigung: " + a)
print("Fehler der aus Messdaten errechneten Beschleunigung: " + da)