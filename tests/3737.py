import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci
import scipy.interpolate as inter

m = np.array([0.76,3.45,7.51,3.5,12.95,3.7,2.2,4.4,-1.46,8.7,0.37,10.7,11.2,2.1,11.4,13.5,11.0,0.0,1.33,3.55,5.2,6,2.2,1,1.35,1,2.6,2.2])
d = np.array([5.1,5.9,5.9,3.6,2.7,3.3,250,4.9,2.7,2.7,3.5,3.5,6.1,75,5.0,2.3,1.3,1.3,1.3,5.7,3.4,3.4,32,130,26,80,190,700])
T = np.array([8100,6030,3480,5120,2700,4600,6250,4700,9600,27000,6580,7000,2950,4430,15000,2400,4130,5800,4130,5400,4100,3850,12500,3400,13500,25000,30000,40000])

[float(i) for i in m]
[float(i) for i in d]
[float(i) for i in T]

M = m + 5 - 5*np.log10(d)

Msun = 4.5
Lsun = 4*(10**33)

Mags = np.array([1.38,11.54,-4.57])
Temps = np.array([9600,27000,3400])

Mscale = Mags/Msun

L = Lsun*(10**((Msun-Mags)/2.5))
sig = 5.67*(10**(-5))

rsq = L/(4*np.pi*sig*(Temps**4))

r = rsq**(0.5)

print '%.2e' %r[0]
print '%.2e' %r[1]
print '%.2e' %r[2]


