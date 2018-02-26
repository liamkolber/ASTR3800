import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick #this library was used to fix the axis problem on the grayscale image
import scipy.stats as sci

#import .fits data as well as the .csv data provided
fitsData = fits.getdata('frame-g-000094-1-0131.fits')
fitsData[fitsData<0.02] = 0.02 #this ignores bad, bright pixel data and makes it dark in order to make actual stars more apparent (value was chosen through guess and check)
fitsData = np.flipud(np.log(fitsData)) #take log of data so it appears on a scale familiar to our eyes, and flipping orients the data correctly
csvData = np.genfromtxt('frame-g-000094-1-0131.csv', delimiter = ",", skip_header = 1) #pull .csv data but skip first line that is just a header
csvData = np.flipud(csvData) #orient data correctly

fig = plt.figure(1) #this figure displays the image given by the fits data
P = fig.add_subplot(111) #subplot was need to make changes to x-axis
plt.imshow(fitsData,cmap = 'Greys_r', extent=(354.473-1025*0.396/3600,354.473+1025*0.396/3600,-0.930-745*0.396/3600,-0.930+745*0.396/3600)) #provides proper coordinates on axis
plt.colorbar()
P.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f')) #makes it so numbers appear properly on the x-axis
plt.show()
