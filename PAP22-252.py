# measure version 1.8.11s
from measure import mv,dsto_mv,val,sqrt,pltext,plt,exp,spcurvefit,tbl,expreg,ln
from numpy import loadtxt,arange,diag

unterg = loadtxt('data/252_untergrund.dat',usecols=[1])

# Silber
n1 = loadtxt('data/252_n1.dat',usecols=[1])
n2 = loadtxt('data/252_n2.dat',usecols=[1])
n3 = loadtxt('data/252_n3.dat',usecols=[1])
n4 = loadtxt('data/252_n4.dat',usecols=[1])

N = n1 + n2 + n3 + n4
N_dsto = sqrt(N)
t = arange(5,405,10)

unterg_ag_mv = mv(4 * unterg)
unterg_ag_mv_dsto = dsto_mv(4*unterg,ddof=0)

def fitfunc(x,A1,l1,A2,l2):
  return A1 * exp(-l1 * x) + A2 * exp(-l2 * x) + unterg_ag_mv

p_opt,p_err = spcurvefit(fitfunc,t,N,p0=[500,0.02,50,0.001],yerr=N_dsto)

Thalb_Ag110 = ln(2) / p_opt[1]
Thalb_Ag110_dsys = ln(2) * p_err[1] / p_opt[1]**2
Thalb_Ag108 = ln(2) / p_opt[3]
Thalb_Ag108_dsys = ln(2) * p_err[3] / p_opt[3]**2

pltext.initplot(num=1,xlabel='Zeit in s',ylabel='# Zerfälle',scale='linlog')
pltext.plotdata(t,N,N_dsto,label='Messwerte',caps=True)
plt.plot(t,fitfunc(t,*p_opt),label='Fit')
pltext.set_layout()

print()
print(val(unterg_ag_mv,unterg_ag_mv_dsto,'Untergrund'))
print(tbl([['','A1','l1','A2','l2'],['Anfangswerte','500','0.02','50','0.001'],['Fitwerte',val(p_opt[0],p_err[0]),val(p_opt[1],p_err[1]),val(p_opt[2],p_err[2]),val(p_opt[3],p_err[3])]]))
print()
print(val(Thalb_Ag108,Thalb_Ag108_dsys,'Ag 108'))
print(val(Thalb_Ag110,Thalb_Ag108_dsys,'Ag 110'))

# Indium
n5 = loadtxt('data/252_n5.dat',usecols=[1])
n5_dsto = sqrt(n5)

t = arange(60,3060,120)

unterg_in_mv = mv(12 * unterg)
unterg_in_mv_dsto = dsto_mv(12 * unterg,ddof=0)

N_in = n5 - unterg_in_mv
N_in_dsto = sqrt(n5_dsto**2 + unterg_in_mv_dsto**2)

pltext.initplot(num=2,scale='linlog')
pltext.plotdata(t,N_in,N_in_dsto,label='Measurements',caps=True)
[slope,dslope] = expreg(t[1:],N_in[1:],N_in_dsto[1:],plot=True)[0:2]
pltext.set_layout()

Thalb_In116 = -ln(2) / slope
Thalb_In116_dsys = ln(2) * dslope / slope**2

print()
print(val(Thalb_In116,Thalb_In116_dsys,'In 116'))

# Plot
print()
plt.show()
