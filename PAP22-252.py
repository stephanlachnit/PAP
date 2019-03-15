# measure version 1.9s
from measure import mv,dsto_mv,val,sqrt,pltext,plt,exp,spcurvefit,tbl,expreg,ln,dev,chi2stat,nplinspace
from numpy import loadtxt,arange

# Untergrund
unterg = loadtxt('data/252_untergrund.dat',usecols=[1])

# Silber
n1 = loadtxt('data/252_n1.dat',usecols=[1])
n2 = loadtxt('data/252_n2.dat',usecols=[1])
n3 = loadtxt('data/252_n3.dat',usecols=[1])
n4 = loadtxt('data/252_n4.dat',usecols=[1])

N_ag = n1 + n2 + n3 + n4
N_ag_err = sqrt(N_ag)
t = arange(5,405,10)

unterg_ag_mv = mv(4 * unterg)
unterg_ag_mv_dsto = dsto_mv(4*unterg,ddof=0)

def fitfunc(x,A1,l1,A2,l2):
  return A1 * exp(-l1 * x) + A2 * exp(-l2 * x) + unterg_ag_mv

def fitfunc_pf(x,A1,l1,A2,l2):
  return A1 * exp(-l1 * x) + A2 * exp(-l2 * x) + unterg_ag_mv + unterg_ag_mv_dsto

def fitfunc_mf(x,A1,l1,A2,l2):
  return A1 * exp(-l1 * x) + A2 * exp(-l2 * x) + unterg_ag_mv - unterg_ag_mv_dsto

p_opt,p_err = spcurvefit(fitfunc,t,N_ag,p0=[500,0.02,50,0.001],yerr=N_ag_err)
p_opt_pf,p_err_pf = spcurvefit(fitfunc_pf,t,N_ag,p0=[500,0.02,50,0.001],yerr=N_ag_err)
p_opt_mf,p_err_mf = spcurvefit(fitfunc_mf,t,N_ag,p0=[500,0.02,50,0.001],yerr=N_ag_err)

Thalb_Ag110 = ln(2) / p_opt[1]
Thalb_Ag110_dsys = ln(2) * p_err[1] / p_opt[1]**2
Thalb_Ag108 = ln(2) / p_opt[3]
Thalb_Ag108_dsys = ln(2) * p_err[3] / p_opt[3]**2

Thalb_Ag110_lit = 24.6
Thalb_Ag108_lit = 2.41 * 60

chi2_ = chi2stat.chi2(N_ag,N_ag_err,fitfunc(t,*p_opt))
chi2_red = chi2stat.chi2_red(chi2_,len(N_ag),ddof=4)
prob = chi2stat.fit_prob(chi2_,len(N_ag),ddof=4)

chi2_pf = chi2stat.chi2(N_ag,N_ag_err,fitfunc_pf(t,*p_opt_pf))
chi2_red_pf = chi2stat.chi2_red(chi2_pf,len(N_ag),ddof=4)
prob_pf = chi2stat.fit_prob(chi2_pf,len(N_ag),ddof=4)

chi2_mf = chi2stat.chi2(N_ag,N_ag_err,fitfunc_mf(t,*p_opt_mf))
chi2_red_mf = chi2stat.chi2_red(chi2_mf,len(N_ag),ddof=4)
prob_mf = chi2stat.fit_prob(chi2_mf,len(N_ag),ddof=4)

t_array = nplinspace(0,400)
pltext.initplot(num=1,title='Abbildung   : Zerfall von Silber',xlabel='Zeit in s',ylabel='# Zerfälle (mit Untergrund)',scale='linlog')
pltext.plotdata(t,N_ag,N_ag_err,label='Messwerte')
plt.plot(t_array,fitfunc(t_array,*p_opt),label='Fit')
plt.plot(t_array,fitfunc_pf(t_array,*p_opt_pf),label='Fit + Fehler Ug')
plt.plot(t_array,fitfunc_mf(t_array,*p_opt_mf),label='Fit - Fehler Ug')
pltext.set_layout(xlim=(0,4e2),ylim=(2e1,4e2))

print('\nSilber:\n')
print(val(unterg_ag_mv,unterg_ag_mv_dsto,name='Untergrund'))
print()
print(tbl([['','A1','l1','A2','l2'],['Fitwerte',val(p_opt[0],p_err[0]),val(p_opt[1],p_err[1]),val(p_opt[2],p_err[2]),val(p_opt[3],p_err[3])],['Fitwerte pf',val(p_opt_pf[0],p_err_pf[0]),val(p_opt_pf[1],p_err_pf[1]),val(p_opt_pf[2],p_err_pf[2]),val(p_opt_pf[3],p_err_pf[3])],['Fitwerte mf',val(p_opt_mf[0],p_err_mf[0]),val(p_opt_mf[1],p_err_mf[1]),val(p_opt_mf[2],p_err_mf[2]),val(p_opt_mf[3],p_err_mf[3])]]))
print()
print(tbl([['','chi2','chi2_red','fitwkeit'],['',val(chi2_),val(chi2_red),val(prob)],['pf',val(chi2_pf),val(chi2_red_pf),val(prob_pf)],['mf',val(chi2_mf),val(chi2_red_mf),val(prob_mf)]]))
print()
print(tbl([['Isotop','Ag 108','Ag 110'],['Thalb (fit)',val(Thalb_Ag108,Thalb_Ag108_dsys),val(Thalb_Ag110,Thalb_Ag110_dsys)],['Thalb (theo)',val(Thalb_Ag108_lit),val(Thalb_Ag110_lit)],dev([Thalb_Ag108,Thalb_Ag110],[Thalb_Ag108_dsys,Thalb_Ag110_dsys],[Thalb_Ag108_lit,Thalb_Ag110_lit],name='Abw',perc=True)]))

# Indium
n5 = loadtxt('data/252_n5.dat',usecols=[1])
n5_dsto = sqrt(n5)

t = arange(60,3060,120)

unterg_in_mv = mv(12 * unterg)
unterg_in_mv_dsto = dsto_mv(12 * unterg,ddof=0)

N_in = n5 - unterg_in_mv
N_in_err = sqrt(n5_dsto**2 + unterg_in_mv_dsto**2)

pltext.initplot(num=2,title='Abblidung   : Zerfall von Indium',xlabel='Zeit in s',ylabel='# Zerfälle (ohne Untergrund)',scale='linlog')
pltext.plotdata(t,N_in,N_in_err,label='Measurements')
[slope,dslope,yitc,dyitc] = expreg(t[1:],N_in[1:],N_in_err[1:],plot=True,prange=(0,3e3))
pltext.set_layout(xlim=(0,3e3),ylim=(3e2,9e2))

Thalb_In116m = -ln(2) / slope
Thalb_In116m_dsys = ln(2) * dslope / slope**2

Thalb_In116m_lit = 54 * 60

chi2_ = chi2stat.chi2(N_in[1:],N_in_err[1:],yitc*exp(slope*t[1:]))
chi2_red = chi2stat.chi2_red(chi2_,len(N_in),ddof=2)
prob = chi2stat.fit_prob(chi2_,len(N_in),ddof=2)

print('\nInduium:\n')
print(val(unterg_in_mv,unterg_in_mv_dsto,name='Untergrund'))
print()
print(tbl([['','A','l'],['Fitwerte',val(yitc,dyitc),val(-slope,dslope)]]))
print()
print(tbl([['chi2','chi2_red','fitwkeit'],[val(chi2_),val(chi2_red),val(prob)]]))
print()
print(val(Thalb_In116m,Thalb_In116m_dsys,'In 116m'))
print(dev(Thalb_In116m,Thalb_In116m_dsys,Thalb_In116m_lit,name='Abw',perc=True))

# Plot
print()
plt.show()
