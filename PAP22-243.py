# measure version 1.8.9s
from measure import pltext,plt,curve_fit,npfarray,val,sqrt,np,nplinspace,mv,T0,kB,dev,dtot_mv
from scipy import integrate

f,U_aus = np.loadtxt('./data/243.txt',skiprows=1,usecols=(0,1),unpack=True)

D=1e-3
U_ein=0.2
g=U_aus/(U_ein*D)

pltext.initplot(num=1,title='Abbildung   : Frequenzgang (Messwerte)',xlabel='Frequenz in Hz',ylabel='g(f)',scale='loglog')
pltext.plotdata(f,g)
pltext.set_layout(xlim=(1e2,1e6),ylim=(3e0,2e3))

def fitfunc(f,V,W1,W2,n1,n2):
  return V/(np.sqrt(1+1/(f/W1)**(2*n1))*np.sqrt(1+(f/W2)**(2*n2)))

V = 1000
W1 = 1000
W2 = 50000
n1 = 5
n2 = 5
p0 = (V,W1,W2,n1,n2)

popt, pcov = curve_fit(fitfunc,f[15:-43],g[15:-43],p0)

print()
print(val(popt[0],sqrt(pcov[0][0]),'V'))
print(val(popt[1],sqrt(pcov[1][1]),'W1'))
print(val(popt[2],sqrt(pcov[2][2]),'W2'))
print(val(popt[3],sqrt(pcov[3][3]),'n1'))
print(val(popt[4],sqrt(pcov[4][4]),'n2'))

pltext.initplot(num=2,title='Abbildung   : Frequenzgang (Fit)',xlabel='Frequenz in Hz',ylabel='g(f)',scale='loglog')
pltext.plotdata(f[17:-45],g[17:-45],label='Messwerte')
f_array = nplinspace(4e2,1.2e5)
plt.plot(f_array,fitfunc(f_array,*popt),label='Fit')
pltext.set_layout(legend=True,xlim=(4e2,1.2e5),ylim=(1e1,2e3))

def fitfuncsquare(f,V,W1,W2,n1,n2):
  return fitfunc(f,V,W1,W2,n1,n2)**2

B = integrate.quad(fitfuncsquare,f[17],f[-45],args=tuple(popt))[0]
B_dsys = 0.02 * B

print(val(B,B_dsys,name='Integral'))

R = npfarray([5,10,15,20,25,30])*1e3
R_dsys = 0.005 * R
U_aus = npfarray([2.4146,3.1253,3.7027,4.2104,4.6677,5.0843])*1e-3
U_aus_err = npfarray([0.00961,0.0113,0.0135,0.0165,0.0182,0.0190])*1e-3
U_V = 1.3841*1e-3
U_V_err = 0.00524*1e-3
D = U_aus**2 - U_V**2
D_err = 2 * sqrt((U_aus * U_aus_err)**2 + (U_V * U_V_err)**2)

def linear(x,c):
  return c*x

popt,pcov=curve_fit(linear,R,D,sigma=D_err)

T_C = npfarray([22.3,23.6])
T_C_dsys = npfarray([0.5,0.5])
T = mv(T_C) + T0
T_dsys = dtot_mv(T_C + T0,T_C_dsys)

c = popt[0]
c_dsto = sqrt(pcov[0][0])

k = c / (4 * T * B)
k_dsto = c_dsto / (4 * T * B)
k_dsys = c / (4 * T * B) * sqrt((T * B_dsys / B)**2 + (B * T_dsys / T)**2)
k_dsys_woT = c / (4 * T * B) * sqrt((T * B_dsys / B)**2)

k_dtot = sqrt(k_dsto**2 + k_dsys**2)

from scipy.stats import chi2

chisqr = np.sum((linear(R,*popt)-D)**2/D_err**2)
dof = 5
chisqr_red = chisqr / dof
prob = 100 * (1-chi2.cdf(chisqr,dof))

pltext.initplot(num=3,title='Abbildung   : Rauschspannung als Funktion des Widerstandes',xlabel='Widerstand in Ohm',ylabel=r'$(U^{2}_{aus}-U^{2}_V)$ / $V^2$')
pltext.plotdata(R,D,D_err,R_dsys,label='Messwerte')
x_array=nplinspace(0,3.2e4)
plt.plot(x_array,linear(x_array,*popt),label='Fitgerade durch Ursprung')
pltext.set_layout(legend=True,xlim=(0,3.5e4),ylim=(0,3e-5))

print(val(T,T_dsys,'Temp'))
print(val(c,c_dsto,'Steigung c'))
print(val(k,name='k'))
print(val(k_dsto,name='k_dsto'))
print(val(k_dsys,name='k_dsys'))
print(val(k_dtot,name='k_dtot'))
print(val(k_dsys_woT,name='k_dsys_woT'))
print(dev(k,k_dtot,kB,name='Abw',perc=True))
print(val(chisqr,name='chi2'))
print(val(chisqr_red,name='chi2_red'))
print(val(prob,name='prob')+'%')

plt.show()
