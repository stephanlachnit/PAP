import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def comma_to_float(valstr):
    return float(valstr.decode("utf-8").replace(',','.'))

def linear (x , a , b):
    return a*x+b

def label( title ,  xlabel , ylabel , grid , xscale , yscale):
    plt.title(title)
    plt.xlabel(xlabel) 
    plt.ylabel(ylabel)
    plt.grid(grid)
    if xscale == '':
        xscale = 'linear'
    plt.xscale(xscale)
    if yscale == '':
        yscale = 'linear'
    plt.yscale(yscale)
    plt.legend()
   
def fenster( nummer ):
    fig = plt.figure(nummer)
    fig.set_size_inches(11.69, 8.27)

def gaus(p , A ,  mu , sigma  ):
     return A * np.exp(-(p-mu)**2/(2*sigma**2)) 

  
x , y = np.loadtxt('Messung.csv', skiprows=2, delimiter = ','  , usecols=(0,1) , converters = {1:comma_to_float, 2:comma_to_float, 3:comma_to_float}, unpack = True) 
y = y - np.mean(y) 

s = 3e-3
ds = 9e-6
T = 23.5
dT = 0.1
alpha = 50e-3
dalpha = 0.05e-3
vgeschw = 0.1e-3
dm = 5
ddm = 0.2
p0 = 101325
T0 = 273.15
TzuPa = 133.22

I = np.array([ 11134 , 11134 , 11133 , 11133 , 11134 ])
dI = np.array([ 105.5 , 105.5 , 105.5 , 105.5 , 105.5])
durchgang1 = np.array([745 , 670 , 590 , 520 , 440 , 370 , 280 , 160 , 70]) * TzuPa
durchgang2 = np.array([750 , 670 , 590 , 520 , 440 , 350 , 280 , 200 , 115 , 35]) * TzuPa
durchgang3 = np.array([745 , 670 , 580 , 515 , 440 , 350 , 270 , 205 , 130 , 60]) * TzuPa
durchgang4 = np.array([745 , 680 , 590 , 520 , 450 , 380 , 305 , 230 , 170 , 80]) * TzuPa
dablese = 10



lamda =  (2)* s / np.mean(I)
dlamda = np.sqrt((2/np.mean(I) * ds)**2 + ((2*s / np.mean(I)**2) * np.mean(dI))**2)



popt , pcov = curve_fit(linear , durchgang1 , np.linspace(dm * (len(durchgang1)-1) , 0 , len(durchgang1)) )
stei1 = popt[0]
achse1 = popt[1]
dstei1 = np.sqrt(pcov[0][0])
popt , pcov = curve_fit(linear , durchgang2 , np.linspace(dm * (len(durchgang2)-1) , 0 , len(durchgang2)) )
stei2 = popt[0]
achse2 = popt[1]
dstei2 = np.sqrt(pcov[0][0])
popt , pcov = curve_fit(linear , durchgang3 , np.linspace(dm * (len(durchgang3)-1) , 0 , len(durchgang3)) )
stei3 = popt[0]
achse3 = popt[1]
dstei3 = np.sqrt(pcov[0][0])
popt , pcov = curve_fit(linear , durchgang4 , np.linspace(dm * (len(durchgang4)-1) , 0 , len(durchgang4)) )
stei4 = popt[0]
achse4 = popt[1]
dstei4 = np.sqrt(pcov[0][0])
stei = (stei1 + stei2 + stei3 + stei4) /4
dstei = np.sqrt(dstei1**2 + dstei2**2 + dstei3**2 + dstei4**2)

n0s = (lamda * p0 * (T + T0)) / (2 * alpha * T0) * stei + 1
dn0s = np.sqrt(((n0s-1)/lamda * dlamda)**2 + ((n0s-1)/T * dT)**2 + ((n0s-1)/alpha * dalpha)**2 + ((n0s-1)/stei * dstei)**2)
"""
p = np.empty(len(durchgang1)+len(durchgang2)+len(durchgang3)+len(durchgang4)-2)
q = 0
while q < len(durchgang1)-1:
       p[q] = durchgang1[q]-durchgang1[q+1]
       q = q+1
q = 0
while q < len(durchgang2)-1:
       p[q+(len(durchgang1)-1)] = durchgang2[q]-durchgang2[q+1]
       q = q+1
q = 0
while q < len(durchgang3)-1:
       p[q+(len(durchgang1)+len(durchgang2)-3)] = durchgang3[q]-durchgang3[q+1]
       q = q+1
q = 0
while q < len(durchgang4)-1:
       p[q+(len(durchgang1)+len(durchgang2)+len(durchgang3)-4)] = durchgang4[q]-durchgang4[q+1]
       q = q+1
dp = np.std(p) / np.sqrt(len(p))
p = np.mean(p)     

n0 = (lamda * dm * p0 * (T + T0)) / (2 * alpha * p * T0) + 1
dn0 = np.sqrt(((n0-1)/lamda * dlamda)**2 + ((n0-1)/T * dT)**2 + ((n0-1)/alpha * dalpha)**2 + ((n0-1)/p * dp)**2)
"""


a = 1.8e-2   # Wert des ersten und letztem Y des Graphen 
i = 1
j = len(y)-1

while y[i] <= a :
    i = i+1
    
while y[j] <= a :
    j = j-1

x = x[i:j]
y = y[i:j]
x = x - np.mean(x)

ypos = np.zeros(len(y))
t = np.zeros(len(y))

f = 0
while f <= len(y)-1:
    ypos[f] = abs(y[f])
    f = f + 1

f = 6 
m = 0
while f <= len(ypos)-7:
    if ypos[f] != np.max(ypos[f-6 : f+6] , axis = 0):
        t[m] = f 
        m = m + 1
    f = f+1

xred = np.delete( x , t , axis = 0)
yposred = np.delete( ypos , t , axis = 0)

init_vals = [0.26199 , 0.0014 , 0.0154]
popt, pcov = curve_fit(gaus, xred, yposred , p0 = init_vals , maxfev = 200000)

L = lamda**2 / (popt[2]*vgeschw*2*np.sqrt(2*np.log(2)))



print('Lambda =' , lamda , '+-' , dlamda)
print('N0 mit Steigung' , n0s , '+-' , dn0s) 
#print('N0 =' , n0 , '+-' , dn0)
print('Kohärenzlänge' , L)

fenster(1)
plt.plot(durchgang1 , np.linspace(dm * (len(durchgang1)-1) , 0 , len(durchgang1)) , marker = '.' , color = 'r' , linewidth = 0 , label = 'Messung 1')
plt.errorbar( durchgang1 , np.linspace(dm * (len(durchgang1)-1) , 0 , len(durchgang1)) , yerr=np.zeros(len(durchgang1))+ddm, xerr= TzuPa *(np.zeros(len(durchgang1))+dablese), color = 'r' , fmt='o', markersize=3)
plt.plot(durchgang1 , linear(durchgang1 , stei1 , achse1) , color = 'g' , label = 'Fit 1')
plt.plot(durchgang2 , np.linspace(dm * (len(durchgang2)-1) , 0 , len(durchgang2)) , marker = '.' , color = 'y' , linewidth = 0 , label = 'Messung 2')
plt.errorbar( durchgang2 , np.linspace(dm * (len(durchgang2)-1) , 0 , len(durchgang2)) , yerr=np.zeros(len(durchgang2))+ddm, xerr= TzuPa *(np.zeros(len(durchgang2))+dablese), color = 'y' , fmt='o', markersize=3)
plt.plot(durchgang2 , linear(durchgang2 , stei2 , achse2) , color = 'b' , label = 'Fit 2')
plt.plot(durchgang3 , np.linspace(dm * (len(durchgang3)-1) , 0 , len(durchgang3)) , marker = '.' , color = 'g' , linewidth = 0 , label = 'Messung 3')
plt.errorbar( durchgang3 , np.linspace(dm * (len(durchgang3)-1) , 0 , len(durchgang3)) , yerr=np.zeros(len(durchgang3))+ddm, xerr= TzuPa *(np.zeros(len(durchgang3))+dablese),  color = 'g' , fmt='o', markersize=3)
plt.plot(durchgang3 , linear(durchgang3 , stei3 , achse3), color = 'r' , label = 'Fit 3')
plt.plot(durchgang4 , np.linspace(dm * (len(durchgang4)-1) , 0 , len(durchgang4)) , marker = '.' , color = 'b' , linewidth = 0 , label = 'Messung 4')
plt.errorbar( durchgang4 , np.linspace(dm * (len(durchgang4)-1) , 0 , len(durchgang4)) , yerr=np.zeros(len(durchgang4))+ddm, xerr= TzuPa *(np.zeros(len(durchgang4))+dablese), color = 'b' , fmt='o', markersize=3)
plt.plot(durchgang4 , linear(durchgang4 , stei4 , achse4), color = 'y' , label = 'Fit 4')
label('Messung des Brechungsindexes ' , 'Druck [Pa]' , 'Anzahl der Interferenzringe' , True , '' , '' )

fenster(2)
plt.plot(x , y , label = 'Messdaten' , color = 'r')
plt.plot(xred, gaus(xred, *popt), label='Fit' , color = 'g')
label('Messung der Kohärenzlänge' , 'Zeit [s]' , 'Spannung [V]' , True , '' , '' )

fenster(3)
plt.plot(x , ypos , color = 'b' )
plt.plot(xred,yposred , linewidth = 1 , color = 'r')
label('positive Messwerte' , 'Zeit' , 'Spannung' , True , '' , '')
