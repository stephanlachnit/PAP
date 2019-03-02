# measure version 1.8.9s
from measure import plt,pltext,npfarray,linreg,sqrt,val

# Aufgabe 1
U = npfarray([420,445,470,495,520,545,570])
N = npfarray([1887,2330,2337,2359,2407,2374,2310])
N_dsto = sqrt(N)

pltext.initplot(num=1,title='Zählrohrcharakteristik',xlabel='Spannung in V',ylabel='Ereignisse')
linreg(U,N,N_dsto,plot=True,frange=range(1,7))
pltext.set_layout()

U0 = 510

# Aufgabe 2

U = [510,610]
N1min = npfarray([9838,9871])
N1min_dsto = sqrt(N1min)
N3min = npfarray([29505,30144])
N3min_dsto = sqrt(N3min)

anstieg_1min = (N1min[1]-N1min[0])
anstieg_1min_dsto = sqrt(N1min_dsto[1]**2 + N1min_dsto[0]**2)
rel_anstieg_1min = anstieg_1min / N1min[0]
rel_anstieg_1min_dsto = 1/N1min[0] * sqrt(anstieg_1min_dsto**2 + (anstieg_1min * N1min_dsto[0] / N1min[0])**2)

anstieg_3min = (N3min[1]-N3min[0])
anstieg_3min_dsto = sqrt(N3min_dsto[1]**2 + N3min_dsto[0]**2)
rel_anstieg_3min = anstieg_3min / N3min[0]
rel_anstieg_3min_dsto = 1/N3min[0] * sqrt(anstieg_3min_dsto**2 + (anstieg_3min * N3min_dsto[0] / N3min[0])**2)


print()
print(val(anstieg_1min,anstieg_1min_dsto,name='Anstieg (1min)'))
print(val(rel_anstieg_1min,rel_anstieg_1min_dsto,name='rel Anstieg (1min)'))
print(val(anstieg_3min,anstieg_3min_dsto,name='Anstieg (3min)'))
print(val(rel_anstieg_3min,rel_anstieg_3min_dsto,name='rel Anstieg (3min)'))

# Aufgabe 3

# Aufgabe 4

# Plot
plt.show()
