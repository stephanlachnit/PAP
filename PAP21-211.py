# measure version 1.3
from measure import pi,sqrt,val,sig

# measured values
l = [0.15,0.25,0.40]
dl = [2e-3,2e-3,2e-3]
f = [[0.616,0.634,0.617,0.633],[0.616,0.665,0.617,0.664],[0.617,0.738,0.616,0.736]]
df = [[6e-3,7e-3,3e-3,5e-3],[6e-3,7e-3,5e-3,6e-3],[6e-3,6e-3,7e-3,8e-3]]

# weak coupling
wI_t = pi * (f[0][1] + f[0][0])
dwI_t = pi * sqrt(df[0][1]**2 + df[0][0]**2)
wII_t = pi * (f[0][1] - f[0][0])
dwII_t = pi * sqrt(df[0][1]**2 + df[0][0]**2)
wI_e = pi * (f[0][3] + f[0][2])
dwI_e = pi * sqrt(df[0][3]**2 + df[0][2]**2)
wII_e = pi * (f[0][3] - f[0][2])
dwII_e = pi * sqrt(df[0][3]**2 + df[0][2]**2)
k_t = (f[0][1]**2 - f[0][0]**2) / (f[0][1]**2 + f[0][0]**2)
dk_t = (4 * f[0][1] * f[0][0]) / (f[0][1]**2 + f[0][0]**2)**2 * sqrt((f[0][0] * df[0][1])**2 + (f[0][1] * df[0][0])**2)
k_e = (f[0][3]**2 - f[0][2]**2) / (f[0][3]**2 + f[0][2]**2)
dk_e = (4 * f[0][3] * f[0][2]) / (f[0][3]**2 + f[0][2]**2)**2 * sqrt((f[0][2] * df[0][3])**2 + (f[0][3] * df[0][2])**2)

print()
print("weak coupling:")
print(val("wI (theo) ", wI_t, dwI_t))
print(val("wI (exp)  ", wI_e, dwI_e))
print(sig("Deviation ", wI_e, dwI_e, wI_t, dwI_t))
print(val("wII (theo)", wII_t, dwII_t))
print(val("wII (exp) ", wII_e, dwII_e))
print(sig("Deviation ", wII_e, dwII_e, wII_t, dwII_t))
print(val("k (theo)  ", k_t, dk_t))
print(val("k (exp)   ", k_e, dk_e))
print(sig("Deviation ", k_e, k_e, k_t, k_t))

# middle coupling
wI_t = pi * (f[1][1] + f[1][0])
dwI_t = pi * sqrt(df[1][1]**2 + df[1][0]**2)
wII_t = pi * (f[1][1] - f[1][0])
dwII_t = pi * sqrt(df[1][1]**2 + df[1][0]**2)
wI_e = pi * (f[1][3] + f[1][2])
dwI_e = pi * sqrt(df[1][3]**2 + df[1][2]**2)
wII_e = pi * (f[1][3] - f[1][2])
dwII_e = pi * sqrt(df[1][3]**2 + df[1][2]**2)
k_t = (f[1][1]**2 - f[1][0]**2) / (f[1][1]**2 + f[1][0]**2)
dk_t = (4 * f[1][1] * f[1][0]) / (f[1][1]**2 + f[1][0]**2)**2 * sqrt((f[1][0] * df[1][1])**2 + (f[1][1] * df[1][0])**2)
k_e = (f[1][3]**2 - f[1][2]**2) / (f[1][3]**2 + f[1][2]**2)
dk_e = (4 * f[1][3] * f[1][2]) / (f[1][3]**2 + f[1][2]**2)**2 * sqrt((f[1][2] * df[1][3])**2 + (f[1][3] * df[1][2])**2)

print()
print("middle coupling:")
print(val("wI (theo) ", wI_t, dwI_t))
print(val("wI (exp)  ", wI_e, dwI_e))
print(sig("Deviation ", wI_e, dwI_e, wI_t, dwI_t))
print(val("wII (theo)", wII_t, dwII_t))
print(val("wII (exp) ", wII_e, dwII_e))
print(sig("Deviation ", wII_e, dwII_e, wII_t, dwII_t))
print(val("k (theo)  ", k_t, dk_t))
print(val("k (exp)   ", k_e, dk_e))
print(sig("Deviation ", k_e, k_e, k_t, k_t))

# strong coupling
wI_t = pi * (f[2][1] + f[2][0])
dwI_t = pi * sqrt(df[2][1]**2 + df[2][0]**2)
wII_t = pi * (f[2][1] - f[2][0])
dwII_t = pi * sqrt(df[2][1]**2 + df[2][0]**2)
wI_e = pi * (f[2][3] + f[2][2])
dwI_e = pi * sqrt(df[2][3]**2 + df[2][2]**2)
wII_e = pi * (f[2][3] - f[2][2])
dwII_e = pi * sqrt(df[2][3]**2 + df[2][2]**2)
k_t = (f[2][1]**2 - f[2][0]**2) / (f[2][1]**2 + f[2][0]**2)
dk_t = (4 * f[2][1] * f[2][0]) / (f[2][1]**2 + f[2][0]**2)**2 * sqrt((f[2][0] * df[2][1])**2 + (f[2][1] * df[2][0])**2)
k_e = (f[2][3]**2 - f[2][2]**2) / (f[2][3]**2 + f[2][2]**2)
dk_e = (4 * f[2][3] * f[2][2]) / (f[2][3]**2 + f[2][2]**2)**2 * sqrt((f[2][2] * df[2][3])**2 + (f[2][3] * df[2][2])**2)

print()
print("strong coupling:")
print(val("wI (theo) ", wI_t, dwI_t))
print(val("wI (exp)  ", wI_e, dwI_e))
print(sig("Deviation ", wI_e, dwI_e, wI_t, dwI_t))
print(val("wII (theo)", wII_t, dwII_t))
print(val("wII (exp) ", wII_e, dwII_e))
print(sig("Deviation ", wII_e, dwII_e, wII_t, dwII_t))
print(val("k (theo)  ", k_t, dk_t))
print(val("k (exp)   ", k_e, dk_e))
print(sig("Deviation ", k_e, k_e, k_t, k_t))
