# measure version 1.8.12s
from measure import mv,dsto_mv,val,sqrt,pltext,plt,exp,spcurvefit,tbl,expreg,ln,dev
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

Thalb_Ag110_lit = 24.6
Thalb_Ag108_lit = 2.41 * 60

pltext.initplot(num=1,title='Abbildung   : Zerfall von Silber',xlabel='Zeit in s',ylabel='# Zerfälle',scale='linlog')
pltext.plotdata(t,N,N_dsto,label='Messwerte')
plt.plot(t,fitfunc(t,*p_opt),label='Fit')
pltext.set_layout(xlim=(0,4e2),ylim=(2e1,4e2))

print()
print(val(unterg_ag_mv,unterg_ag_mv_dsto,'Untergrund'))
print(tbl([['','A1','l1','A2','l2'],['Anfangswerte','500','0.02','50','0.001'],['Fitwerte',val(p_opt[0],p_err[0]),val(p_opt[1],p_err[1]),val(p_opt[2],p_err[2]),val(p_opt[3],p_err[3])]]))
print()
print(val(Thalb_Ag108,Thalb_Ag108_dsys,'Ag 108'))
print(dev(Thalb_Ag108,Thalb_Ag108_dsys,Thalb_Ag108_lit,name='Abw',perc=True))
print(val(Thalb_Ag110,Thalb_Ag110_dsys,'Ag 110'))
print(dev(Thalb_Ag110,Thalb_Ag110_dsys,Thalb_Ag110_lit,name='Abw',perc=True))

# Indium
n5 = loadtxt('data/252_n5.dat',usecols=[1])
n5_dsto = sqrt(n5)

t = arange(60,3060,120)

unterg_in_mv = mv(12 * unterg)
unterg_in_mv_dsto = dsto_mv(12 * unterg,ddof=0)

N_in = n5 - unterg_in_mv
N_in_dsto = sqrt(n5_dsto**2 + unterg_in_mv_dsto**2)

pltext.initplot(num=2,title='Abblidung   : Zerfall von Indium',xlabel='Zeit in s',ylabel='# Zerfälle',scale='linlog')
pltext.plotdata(t,N_in,N_in_dsto,label='Measurements')
[slope,dslope,yitc,dyitc] = expreg(t[1:],N_in[1:],N_in_dsto[1:],plot=True,prange=(0,3e3))
pltext.set_layout(xlim=(0,3e3),ylim=(3e2,9e2))

Thalb_In116m = -ln(2) / slope
Thalb_In116m_dsys = ln(2) * dslope / slope**2

Thalb_In116m_lit = 54 * 60

print()
print(tbl([['','A','l'],['Fitwerte',val(yitc,dyitc),val(-slope,dslope)]]))
print()
print(val(Thalb_In116m,Thalb_In116m_dsys,'In 116m'))
print(dev(Thalb_In116m,Thalb_In116m_dsys,Thalb_In116m_lit,name='Abw',perc=True))

# Plot
print()
plt.show()
