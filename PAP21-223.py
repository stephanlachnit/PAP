# measure version 1.8.2
from measure import np,plt,pltext,npfarray,curve_fit,mv,dsys_mv,dsto_mv,dtot,dsto,val,sqrt,pi,kB,sig,T0
from scipy.stats import norm

# Messwerte
t = npfarray([2.001,3.001,4.001,5.002,6.002,7.003,8.003,9.105,10.105,11.106,12.106,13.107,14.107,15.208,16.208,17.309,18.309,19.41,20.41,21.511,22.511,23.612,24.612,25.712,26.714,27.814,28.815,29.915,30.916,32.016,33.017,34.117,35.118,36.218,37.219,38.319,39.32,40.42,41.421,42.521,43.522,44.623,45.624,46.724,47.725,48.825,49.826,50.926,51.927,53.027,54.028,55.128,56.129,57.229,58.23,59.33,60.331,61.432,62.433,63.533,64.534,65.634,66.635,67.735,68.736,69.836,70.837,71.937,72.938,74.038,75.039,76.139,77.14,78.24,79.242,80.342,81.343,82.443,83.444,84.544,85.545,86.645,87.646,88.746,89.747,90.847,91.847,92.948,93.948,95.048,96.049,97.151,98.151,99.252,100.252,101.253,102.253,103.254,104.254,105.255,106.255,107.256,108.256,109.256,110.257,111.257,112.258,113.258,114.26,115.26,116.261,117.261,118.262,119.262,120.363,121.363,122.464,123.464,124.565,125.565,126.666,127.666,128.767,129.767,130.868,131.869,132.97,133.97,135.071,136.071,137.172,138.172,139.273,140.273,141.374,142.374,143.475,144.475,145.576,146.576,147.676,148.677,149.779,150.779,151.88,152.88,153.981,154.981,156.081,157.082,158.182,159.183,160.283,161.284,162.384,163.385,164.485,165.486,166.587,167.588,168.688,169.689,170.789,171.79,172.89,173.891,174.991,175.992,177.089])
x = 1e-6 * npfarray([59.849,59.304,59.576,59.122,59.394,58.667,58.213,57.667,58.122,58.758,58.849,57.394,57.122,56.485,57.485,57.576,57.667,57.304,55.849,56.031,56.031,56.758,56.485,55.667,55.576,55.667,55.485,56.213,57.394,57.667,58.758,60.122,60.304,59.031,58.304,57.031,56.122,54.304,54.849,55.304,55.485,55.758,55.122,54.485,55.304,55.304,55.758,57.213,57.576,54.758,53.849,53.576,53.758,53.213,53.394,53.394,53.122,52.849,52.031,51.304,52.304,53.667,53.213,52.667,53.758,53.031,51.758,51.394,51.485,52.031,53.304,53.304,53.667,52.849,53.304,52.576,52.94,51.122,51.758,51.213,50.576,51.758,51.94,51.213,51.122,48.667,47.485,47.849,47.122,46.394,45.94,46.94,46.576,46.94,47.213,46.394,46.213,46.576,46.485,45.213,45.485,45.213,45.031,43.94,44.394,43.394,43.667,44.667,44.304,43.394,43.213,42.304,44.304,44.849,44.213,44.213,44.394,44.576,45.122,44.213,42.94,43.485,43.758,44.485,43.667,44.122,44.849,44.667,43.94,43.758,44.122,43.758,42.758,43.485,44.031,43.94,43.758,42.849,44.213,44.576,44.122,43.849,45.031,45.394,44.758,45.667,46.304,46.667,46.667,47.304,47.485,47.304,46.213,46.849,46.94,45.576,45.758,45.849,46.94,46.213,45.213,45.394,46.304,46.122,47.849,47.94,48.394,49.94,49.213])
y = 1e-6 * npfarray([50.821,52.548,53.185,52.094,52.276,53.639,54.548,56.003,56.457,55.821,56.003,56.185,55.548,56.366,57.003,57.003,57.912,58.094,58.457,58.548,57.639,57.548,57.457,56.366,56.457,55.821,56.548,56.639,56.548,56.821,56.094,56.094,56.639,56.912,56.366,56.276,57.094,57.457,58.548,59.276,59.639,61.639,62.639,63.457,64.639,63.912,63.457,63.639,63.548,63.912,63.548,65.73,66.366,67.094,66.73,66.73,67.094,67.276,68.276,68.185,69.094,68.548,67.912,67.276,67.366,66.821,66.821,65.73,65.457,66.457,67.003,67.548,68.821,68.821,70.094,67.821,69.821,70.639,70.094,70.73,70.457,69.821,69.73,69.821,69.457,69.185,69.185,69.185,70.185,70.276,69.912,70.185,71.821,73.094,73.821,72.639,72.639,73.548,73.366,72.912,73.548,72.821,70.912,70.185,70.912,71.639,70.912,69.821,71.276,70.548,71.185,71.366,71.548,71.548,70.457,69.094,69.276,70.457,71.003,72.185,72.548,72.548,73.548,73.548,74.548,74.548,74.366,76.094,76.73,76.185,76.094,75.276,75.73,75.548,75.548,75.457,75.457,75.366,77.548,77.821,79.094,79.73,79.639,80.003,79.821,79.457,79.912,78.912,77.457,76.094,75.366,75.821,76.003,74.821,75.366,74.457,74.185,74.73,76.003,75.003,75.003,74.003,74.73,75.276,75.639,75.912,75.73,76.185,77.548])

dt = dx = dy = npfarray([])
for n in range(168):
  dt = np.append(dt, t[n+1] - t[n])
  dx = np.append(dx, x[n+1] - x[n])
  dy = np.append(dy, y[n+1] - y[n])

r_k = 755e-9 / 2.
r_k_dsys = 30e-9 / 2.
T = npfarray([22.6,23.0])
T_dsys = npfarray([0.1,0.1])
T_mv = mv(T) + T0
T_dtot = dtot(dsys_mv(T_dsys),dsto_mv(T))
nu = 9.40e-4
nu_dsys = 0.05e-4

# Teilchenbewegung
pltext.initplot(num=1,title='Bewegung des Teilchens',xlabel='x in m',ylabel='y in m')
plt.plot(x,y, marker='s')
plt.savefig('fig1.pdf', format='pdf')

# Histogramm
r_sqr = (dx**2+dy**2)
r_sqr_mv = mv(r_sqr)
r_sqr_dsto = dsto_mv(r_sqr)
dt_mv = mv(dt)
dt_dsto = dsto_mv(dt)

hist_D = r_sqr_mv / (4. * dt_mv)
hist_D_dtot = 1. / (4. * dt_mv) * sqrt(r_sqr_dsto**2 + (r_sqr_mv * dt_dsto / dt_mv)**2)
hist_kB = 6.*pi * nu * r_k * hist_D / T_mv
hist_kB_dtot = 6.*pi / T_mv * sqrt((nu_dsys * r_k * hist_D)**2 + (nu * r_k_dsys * hist_D)**2 + (nu * r_k * hist_D_dtot)**2 + (nu * r_k * hist_D * T_dtot / T_mv)**2)

d_all = 1e6 * np.append(dx,dy)
points = npfarray([0.01 * n for n in range(-300,301)])
mu = mv(d_all)
sigma = dsto(d_all)
gauss = norm.pdf(points, mu, sigma)

pltext.initplot(num=2,title='Histogramm der Verschiebungen',xlabel='Verschiebung in μm',ylabel='relative Häufigkeit')
plt.xlim(-3,3)
plt.ylim(0,0.6)
plt.hist(d_all,density=True)
plt.plot(points,gauss)
plt.text(-2.41,0.504,'σ = '+val('',sigma)+'\nμ = '+val('',mu))
plt.savefig('fig2.pdf', format='pdf')

print()
print('Histogramm:')
print(val('Diffusionskonstante D', hist_D, hist_D_dtot))
print(val('Bolatzmannkonstante k', hist_kB, hist_kB_dtot))

# Kumulative Verteilung der Verschiebungsquadrate
r_kum = np.cumsum(r_sqr)
def lin(x,a,b):
  return a*x+b
popt,pcov = curve_fit(lin, t[:-1], r_kum)
slope = popt[0]
yitc = popt[1]
slope_dtot = sqrt(pcov[0][0])
yitc_dtot = sqrt(pcov[1][1])

vert_D = slope / 4.
vert_D_dtot = slope_dtot / 4.
vert_kB = 6.*pi * nu * r_k * vert_D / T_mv
vert_kB_dtot = 6.*pi / T_mv * sqrt((nu_dsys * r_k * vert_D)**2 + (nu * r_k_dsys * vert_D)**2 + (nu * r_k * vert_D_dtot)**2 + (nu * r_k * vert_D * T_dtot / T_mv)**2)

pltext.initplot(num=3,title='Kumulative Verteilung der Verschiebungsquadrate',xlabel='Zeit in s',ylabel='Summe $r_i^2$ in μm')
pltext.plotdata(t[:-1], r_kum)
plt.plot(t[:-1], lin(t[:-1],slope,yitc), label='Ausgleichsgerade')
plt.plot(t[:-1], lin(t[:-1],slope-slope_dtot,yitc+yitc_dtot), label='Fehlergerade')
plt.xlim(0.0,177.0)
plt.ylim(-0.04e-10,2.20e-10)
plt.legend(loc='upper left')
plt.savefig('fig3.pdf', format='pdf')

print()
print('Kumulative Verteilung:')
print(val('Diffusionskonstante D', vert_D, vert_D_dtot))
print(val('Bolatzmannkonstante k', vert_kB, vert_kB_dtot))

# Vergleich der Werte
print()
print('Vergleich der Werte zueinander:')
print(sig('Abweichung',hist_kB,hist_kB_dtot,vert_kB,vert_kB_dtot))
print()
print('Vergleich mit dem Literaturwert:')
print(sig('Histogramm',hist_kB,hist_kB_dtot,kB))
print(sig('Verteilung',vert_kB,vert_kB_dtot,kB))

plt.show()
