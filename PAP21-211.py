# measure version 1.3
from measure import pi,mean_value,std_dev_m,sqrt,val,sig,linreg,plot

# measured values
l = [0.15,0.25,0.40]
dl = [2e-3,2e-3,2e-3]
f = [[0.616,0.634,0.617,0.633],[0.616,0.665,0.617,0.664],[0.617,0.738,0.616,0.736]]
df = [[6e-3,7e-3,3e-3,5e-3],[6e-3,7e-3,5e-3,6e-3],[6e-3,6e-3,7e-3,8e-3]]
T = [1.62, 1.63, 1.62,1.63,1.63,1.63]

# no coupling
w = 2 * pi / mean_value(T)
dw = 2 * pi * std_dev_m(T) / mean_value(T)**2

print()
print("no coupling;")
print(val("w",w, dw))

# weak coupling
w1_w = 2 * pi * f[0][0]
dw1_w = 2 * pi * df[0][0]
w2_w = 2 * pi * f[0][1]
dw2_w = 2 * pi * df[0][1]
wI_t_w = pi * (f[0][1] + f[0][0])
dwI_t_w = pi * sqrt(df[0][1]**2 + df[0][0]**2)
wII_t_w = pi * (f[0][1] - f[0][0])
dwII_t_w = pi * sqrt(df[0][1]**2 + df[0][0]**2)
wI_e_w = pi * (f[0][3] + f[0][2])
dwI_e_w = pi * sqrt(df[0][3]**2 + df[0][2]**2)
wII_e_w = pi * (f[0][3] - f[0][2])
dwII_e_w = pi * sqrt(df[0][3]**2 + df[0][2]**2)
k_t_w = (f[0][1]**2 - f[0][0]**2) / (f[0][1]**2 + f[0][0]**2)
dk_t_w = (4 * f[0][1] * f[0][0]) / (f[0][1]**2 + f[0][0]**2)**2 * sqrt((f[0][0] * df[0][1])**2 + (f[0][1] * df[0][0])**2)
k_e_w = (f[0][3]**2 - f[0][2]**2) / (f[0][3]**2 + f[0][2]**2)
dk_e_w = (4 * f[0][3] * f[0][2]) / (f[0][3]**2 + f[0][2]**2)**2 * sqrt((f[0][2] * df[0][3])**2 + (f[0][3] * df[0][2])**2)

print()
print("weak coupling:")
print(val("w1        ", w1_w, dw1_w))
print(val("w2        ", w2_w, dw2_w))
print(val("wI (theo) ", wI_t_w, dwI_t_w))
print(val("wI (exp)  ", wI_e_w, dwI_e_w))
print(sig("Deviation ", wI_e_w, dwI_e_w, wI_t_w, dwI_t_w))
print(val("wII (theo)", wII_t_w, dwII_t_w))
print(val("wII (exp) ", wII_e_w, dwII_e_w))
print(sig("Deviation ", wII_e_w, dwII_e_w, wII_t_w, dwII_t_w))
print(val("k (theo)  ", k_t_w, dk_t_w))
print(val("k (exp)   ", k_e_w, dk_e_w))
print(sig("Deviation ", k_e_w, dk_e_w, k_t_w, dk_t_w))

# middle coupling
w1_m = 2 * pi * f[1][0]
dw1_m = 2 * pi * df[1][0]
w2_m = 2 * pi * f[1][1]
dw2_m = 2 * pi * df[1][1]
wI_t_m = pi * (f[1][1] + f[1][0])
dwI_t_m = pi * sqrt(df[1][1]**2 + df[1][0]**2)
wII_t_m = pi * (f[1][1] - f[1][0])
dwII_t_m = pi * sqrt(df[1][1]**2 + df[1][0]**2)
wI_e_m = pi * (f[1][3] + f[1][2])
dwI_e_m = pi * sqrt(df[1][3]**2 + df[1][2]**2)
wII_e_m = pi * (f[1][3] - f[1][2])
dwII_e_m = pi * sqrt(df[1][3]**2 + df[1][2]**2)
k_t_m = (f[1][1]**2 - f[1][0]**2) / (f[1][1]**2 + f[1][0]**2)
dk_t_m = (4 * f[1][1] * f[1][0]) / (f[1][1]**2 + f[1][0]**2)**2 * sqrt((f[1][0] * df[1][1])**2 + (f[1][1] * df[1][0])**2)
k_e_m = (f[1][3]**2 - f[1][2]**2) / (f[1][3]**2 + f[1][2]**2)
dk_e_m = (4 * f[1][3] * f[1][2]) / (f[1][3]**2 + f[1][2]**2)**2 * sqrt((f[1][2] * df[1][3])**2 + (f[1][3] * df[1][2])**2)

print()
print("middle coupling:")
print(val("w1        ", w1_m, dw1_m))
print(val("w2        ", w2_m, dw2_m))
print(val("wI (theo) ", wI_t_m, dwI_t_m))
print(val("wI (exp)  ", wI_e_m, dwI_e_m))
print(sig("Deviation ", wI_e_m, dwI_e_m, wI_t_m, dwI_t_m))
print(val("wII (theo)", wII_t_m, dwII_t_m))
print(val("wII (exp) ", wII_e_m, dwII_e_m))
print(sig("Deviation ", wII_e_m, dwII_e_m, wII_t_m, dwII_t_m))
print(val("k (theo)  ", k_t_m, dk_t_m))
print(val("k (exp)   ", k_e_m, dk_e_m))
print(sig("Deviation ", k_e_m, dk_e_m, k_t_m, dk_t_m))

# strong coupling
w1_s = 2 * pi * f[2][0]
dw1_s = 2 * pi * df[2][0]
w2_s = 2 * pi * f[2][1]
dw2_s = 2 * pi * df[2][1]
wI_t_s = pi * (f[2][1] + f[2][0])
dwI_t_s = pi * sqrt(df[2][1]**2 + df[2][0]**2)
wII_t_s = pi * (f[2][1] - f[2][0])
dwII_t_s = pi * sqrt(df[2][1]**2 + df[2][0]**2)
wI_e_s = pi * (f[2][3] + f[2][2])
dwI_e_s = pi * sqrt(df[2][3]**2 + df[2][2]**2)
wII_e_s = pi * (f[2][3] - f[2][2])
dwII_e_s = pi * sqrt(df[2][3]**2 + df[2][2]**2)
k_t_s = (f[2][1]**2 - f[2][0]**2) / (f[2][1]**2 + f[2][0]**2)
dk_t_s = (4 * f[2][1] * f[2][0]) / (f[2][1]**2 + f[2][0]**2)**2 * sqrt((f[2][0] * df[2][1])**2 + (f[2][1] * df[2][0])**2)
k_e_s = (f[2][3]**2 - f[2][2]**2) / (f[2][3]**2 + f[2][2]**2)
dk_e_s = (4 * f[2][3] * f[2][2]) / (f[2][3]**2 + f[2][2]**2)**2 * sqrt((f[2][2] * df[2][3])**2 + (f[2][3] * df[2][2])**2)

print()
print("strong coupling:")
print(val("w1        ", w1_s, dw1_s))
print(val("w2        ", w2_s, dw2_s))
print(val("wI (theo) ", wI_t_s, dwI_t_s))
print(val("wI (exp)  ", wI_e_s, dwI_e_s))
print(sig("Deviation ", wI_e_s, dwI_e_s, wI_t_s, dwI_t_s))
print(val("wII (theo)", wII_t_s, dwII_t_s))
print(val("wII (exp) ", wII_e_s, dwII_e_s))
print(sig("Deviation ", wII_e_s, dwII_e_s, wII_t_s, dwII_t_s))
print(val("k (theo)  ", k_t_s, dk_t_s))
print(val("k (exp)   ", k_e_s, dk_e_s))
print(sig("Deviation ", k_e_s, dk_e_s, k_t_s, dk_t_s))

# proportionality between coupling and lengths
l2 = [l[i]**2 for i in range(len(l))]
dl2 = [2 * l[i] * dl[i] for i in range(len(l))]
k = [k_t_w, k_t_m, k_t_s]
dk = [dk_t_w, dk_t_m, dk_t_s]
[g, dg, b, db, plt] = linreg(l2, k, dk, dl2, drawplot=True, xlabel="l^2 / m^2", ylabel="coupling constant k")
plt.saveplot("linreg.pdf")

print()
print("proportionality:")
print(val("slope", g, dg))
print(val("y-itc", b, db))

# ratios between k's and l^2's
rk_wm = k_e_w / k_e_m
drk_wm = 1 / k_e_m * sqrt(dk_e_w**2 + (k_e_w / k_e_m * dk_e_m)**2)
rk_ws = k_e_w / k_e_s
drk_ws = 1 / k_e_s * sqrt(dk_e_w**2 + (k_e_w / k_e_s * dk_e_s)**2)
rk_ms = k_e_m / k_e_s
drk_ms = 1 / k_e_s * sqrt(dk_e_m**2 + (k_e_m / k_e_s * dk_e_s)**2)
l2_wm = l2[0] / l2[1]
dl2_wm = 1 / l2[1] * sqrt(dl2[0]**2 + (l2[0] / l2[1] * dl2[1])**2)
l2_ws = l2[0] / l2[2]
dl2_ws = 1 / l2[2] * sqrt(dl2[0]**2 + (l2[0] / l2[2] * dl2[2])**2)
l2_ms = l2[1] / l2[2]
dl2_ms = 1 / l2[2] * sqrt(dl2[1]**2 + (l2[1] / l2[2] * dl2[2])**2)

print(val("kw / km  ", rk_wm, drk_wm))
print(val("l2w / l2m", l2_wm, dl2_wm))
print(sig("deviation", rk_wm, drk_wm, l2_wm, dl2_wm))
print(val("kw / ks  ", rk_ws, drk_ws))
print(val("l2w / l2s", l2_ws, dl2_ws))
print(sig("deviation", rk_ws, drk_ws, l2_ws, dl2_ws))
print(val("km / ks  ", rk_ms, drk_ms))
print(val("l2m / l2s", l2_ms, dl2_ms))
print(sig("deviation", rk_ms, drk_ms, l2_ms, dl2_ms))

plot.showplots()
