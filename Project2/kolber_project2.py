import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick #this library was used to fix the axis problem on the grayscale image
import scipy.stats as sci

def basicPlot(ttl,x,y,s): #streamline plot labeling
	plt.title(ttl,fontsize=16)
	plt.xlabel(x,fontsize=14)
	plt.ylabel(y,fontsize=14)
	if s == True: #True if user wants to show plot
		plt.show() 
#-----------------------------------------------------------------------
print "\nBeginning Project"
print "Press enter to begin Problem A(i)"
my.pause()
#-----------------------------------------------------------------------

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
basicPlot('Grayscale Image',x = 'Right Ascension [degrees]',y = 'Declination [degrees]',s = True)

#-----------------------------------------------------------------------
print "End of Problem A(i)"
print "Press enter to begin Problem A(ii) and A(iii)"
my.pause()
#-----------------------------------------------------------------------
#this section splits .csv data and plot

picX = []
picY = []
picG = []
for i in csvData: #this loop splits data
	picY.append(i[0])
	picX.append(i[1])
	temp = -(i[2]**3)/100+200 #edit data so that it properly shows sizes on scatter plot (odd scaling is result of guess and check)
	picG.append(temp) #apply scale to g-band so that differences in magnitude are more obvious

plt.figure(2) #this figure plots a scatter plot of all the stars with scaled g-band
plt.scatter(picX, picY, s = picG, marker = "*", alpha = 0.5)
plt.xlim(0,2048)
plt.ylim(0,1498)
basicPlot('Scatter Plot of Star Locations (Size Reflects Brightness)',x = 'Right Ascension [pixels]',y = 'Declination [pixels]',s = True)

#-----------------------------------------------------------------------
print "End of Problem A(ii) and A(iii)"
print "Press enter to begin Problem B(i)"
my.pause()
#-----------------------------------------------------------------------
#this section generates random data and plot in a scatter plot

amt = len(csvData)
seed = np.random.seed(414928)
genX = np.random.uniform(0,2048,amt)
genY = np.random.uniform(0,1498,amt)
randData = []
for i in range(amt): #this loop ensures that data is properly formatted in list for easy access later
	randData.append([genX[i],genY[i]])

plt.figure(3) #this figure plots a scatter plot of the randomly generated data
plt.scatter(genX, genY, s = 40, marker = "*", alpha = 0.5)
plt.xlim(0,2048)
plt.ylim(0,1498)
basicPlot('Scatter Plot of Star Locations (Generated Randomly)',x = 'Right Ascension [pixels]',y = 'Declination [pixels]',s = True)

#-----------------------------------------------------------------------
print "End of Problem B(i)"
print "Press enter to begin Problem B(ii,a)"
my.pause()
#-----------------------------------------------------------------------
#this section calculates each stars distance to their closest neightbor

def findMin(data,x,y,limit): #calculate minimum values
	minVal = np.zeros(limit)
	index = 0
	for star in data:
		dx = x - star[0]
		dy = y - star[1]
		dist = np.sort(np.sqrt(dx**2 + dy**2)) #sort data after distance calculation so that minimum is first in list
		minVal[index] = dist[1] #take second value because first value is zero (star's distance from itself)
		index += 1
	return minVal

minD = findMin(csvData,x=picY,y=picX,limit=amt)
randminD = findMin(randData,x=genX,y=genY,limit=amt)

plt.figure(4) #this figure plots histograms of minimum data for comparison
plt.subplot(211)
plt.hist(minD,10,range=(0,randminD.max())) #histogram of .csv data
basicPlot("Distribution of Minimum Distances for Given (top plot) and Randomly Generated (bottom plot) Data",x = "Bin Edges",y = "Amount per Bin",s = False)
plt.figure(4)
plt.subplot(212)
plt.hist(randminD,10,range=(0,randminD.max())) #histogram of random data
basicPlot("",x = "Bin Edges",y = "Amount per Bin",s = True)

#-----------------------------------------------------------------------
print "End of Problem B(ii,a)"
print "Press enter to begin Problem B(ii,b)"
my.pause()
#-----------------------------------------------------------------------
#this section calculates a X_squared value and probability to exceed

def Xsq(binEx,binOb): #calculate X_squared(expected,observed)
	tot = 0
	temp = []
	
	for i in range(len(binOb)-1):
		if binOb[i] >= 5 and binEx[i] >= 5: #this checks that there are at least 5 values per bin
			tot = ((binOb[i] - (binEx[i]*1.0))**2)/binEx[i]
			temp.append(tot)
	Sum = np.sum(temp)
	
	return Sum

def pval(nbins,chi): #calculate probability to exceed (number of bins, chi_squared)
	p = 1 - sci.chi2.cdf(x = chi,df = nbins - 1)
	return p

#get bin sizes for expected and observed data
expected = np.histogram(randminD)[0] #random data
observed = np.histogram(minD)[0] #.csv data

print "X_squared of picture data: ",Xsq(expected,observed)
print "Probability to exceed of picture data: ",pval(10,Xsq(expected,observed)),"\n"

#-----------------------------------------------------------------------
print "End of Problem B(ii,b)"
print "Press enter to begin Problem B(ii,c)"
my.pause()
#-----------------------------------------------------------------------
#this section calculates X_squared and probability to exceed ignoring the closest neighbors (the first bin)

print "X_squared ignoring closest neightbors: ",Xsq(expected[1:],observed[1:])
print "Probability to exceed ignoring closest neighbors: ",pval(9,Xsq(expected[1:],observed[1:])),"\n"

#-----------------------------------------------------------------------
print "End of Problem B(ii,c)"
print "Press enter to end project"
my.pause()
#-----------------------------------------------------------------------



