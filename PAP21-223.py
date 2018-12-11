import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.optimize import curve_fit

def comma_to_float(valstr):
    return float(valstr.decode("utf-8").replace(',','.'))

def linear (x , a , b):
    return a*x+b

def label( title ,  xlabel , ylabel , grid , xscale , yscale ):
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
   
def fenster( nummer ):
    fig = plt.figure(nummer)
    plt.clf()
    fig.set_size_inches(11.69, 8.27)
    
t,x,y = np.loadtxt('data/Messung.dat', skiprows=1, usecols=(1,2,3), converters = {1:comma_to_float, 2:comma_to_float, 3:comma_to_float}, unpack = True)

d_kugel = 755e-9
dd_kugel = 30e-9
nü = 9.5e-4
dnü = 0.07e-4
Temp = 22.8+273.15
dTemp = np.sqrt(0.1**2 + 0.1**2)
dt = np.array([])
dx = np.array([])
dy= np.array([])
i=0
while i < len(t)-1:
    dt = np.append(dt, t[i+1]-t[i])
    dx = np.append(dx, x[i+1]-x[i])
    dy = np.append(dy, y[i+1]-y[i])
    i = i + 1

dx *= 1e-6
dy *= 1e-6
r_squared = (dx**2+dy**2)
r_squared_mean = np.mean(r_squared)
r_squared_mean_std = np.std(r_squared) / np.sqrt(len(r_squared))
dt_mean = np.mean(dt)
all_data = np.append(dx , dy)
mu = np.mean(all_data)
sigma = np.std(all_data)
gauss = mlab.normpdf(np.linspace(-4,4,20000), mu , sigma)
r_kumm = np.cumsum(r_squared)
popt , pcov = curve_fit(linear , t[:-1] , r_kumm)
steigung = popt[0]
dsteigung = np.sqrt(pcov[0][0]) #?
k_bolz = (6*np.pi*nü*(d_kugel/2)*r_squared_mean)/(4*Temp*dt_mean)
D = r_squared_mean  / (4*dt_mean)
k_bolz2 = (6*np.pi*nü*(d_kugel/2)*steigung)/(4*Temp)
D2 = (steigung)/4

print("r_squared_mean=" , r_squared_mean)
print("r_squared_mean_std =" , r_squared_mean_std)
print("dt_mean= " , dt_mean)
print('steigung',steigung,' +- ',dsteigung)
print("Bolzmannkonstante " , k_bolz)
print("Diffusionskonstante " , D)
print("Bolzmannkonstante mit Steigung " , k_bolz2)
print("Diffusionskonstante mit Steigung " , D2)


fenster( 1 )
plt.plot(x,y, marker = 's', color = 'red', linewidth = 1)
label('Brownsche Bewegung' , 'x ['+r'$\mu$'+'m'']', 'y ['+r'$\mu$'+'m'']' , False , '' , '' )
plt.savefig('figures/Graph1.pdf' , format = 'PDF')

fenster( 2 )
plt.hist(all_data , normed = 1)
plt.plot(np.linspace(-4 , 4 , 20000) , gauss , 'b-' , linewidth=2)
label('Histogramm der Verschiebungen' , 'Verschiebung [$\\mu m$]' , 'rel. Häufigkeit' , True , '' , '')
plt.savefig('figures/Graph2.pdf' , format = 'PDF')

fenster( 3 )
plt.plot(t[:-1] , r_kumm , marker='.' , color='red' , linewidth=0)
plt.plot(t[:-1] , linear(t[:-1] , *popt))
label('Kummulative Verschiebung' , 'Zeit [s]' , 'Summe $r_i^2 [\\mu m^2$]' , True , '' , '')
plt.savefig('figures/Graph3.pdf' , format = 'PDF')
