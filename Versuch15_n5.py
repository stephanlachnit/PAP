import measure.py

g = 

m = 

h1 = 
h2 = 

t1 = []
t2 = []
s1 = []
s2 = []

t1m = mean_value(t1)
t2m = mean_value(t2)
s1m = mean_value(s1)
s2m = mean_value(s2)

ek = 0.5 * m * ((s2m - s1m) / (t2m - t1m))**2
dek = 0.5 * m * sqrt((2.0 * (s2m - s1m)**2 / (t2m - t1m)**3 * std_dev_m(t1))**2
                     + (2.0 * (s2m - s1m)**2 / (t2m - t1m)**3 * std_dev_m(t2))**2
                     + (2.0 * (s2m - s1m) / (t2m - t1m)**2 * std_dev_m(s1))**2
                     + (2.0 * (s2m - s1m) / (t2m - t1m)**2 * std_dev_m(s2))**2)

ep = m * g