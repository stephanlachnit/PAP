import measure.py

t_m = [[], [], [], []]
s = [ , , , ]
ds = [ , , , ]

t = [mean_value(x) for x in t_m]
dt = [std_dev_m(x) for x in t_m]

t2 = [x**2 for x in t]
dt2 = [abs(2*t[i] * dt[i]) for i in range(len(dt))]


