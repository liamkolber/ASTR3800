import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci

my.introduction('21 Sep, 2016','30 Sep, 2016')
#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin Problem A(i)"
my.pause()
#-----------------------------------------------------------------------
#this section reads in file data and plots both Red and Ca II rdc and contrast images

def readFitsIn(time,color,typ): #this function reads in the fits data
	data = fits.getdata("20050702."+str(time)+".HW."+color+".P."+typ+".fits")
	return data

def specTime(time): #this function accesses all the fits data available for the given time
	Rrdc = readFitsIn(time,'R','rdc')
	Rcont = readFitsIn(time,'R','contrast')
	Krdc = readFitsIn(time,'K','rdc')
	Kcont = readFitsIn(time,'K','contrast')
	return [Rrdc,Rcont,Krdc,Kcont]

'''
David: These functions provide a very compact way of handling the file names, which is good.
The downside of doing it this way is that every time you need access to a FITS file, you have
to read in every other FITS file of the same observation time. This creates a lot of unnecessary
I/O.
'''

plt.figure(1) #this figure contains the four solar images
for i in range(len(specTime(1702))):
	sp = plt.subplot(221+i)
	plt.imshow(specTime(1702)[i],cmap='Greys_r')
	sp.set_xticklabels([]) #removes axis labels
	sp.set_yticklabels([]) #removes axis labels
plt.suptitle('Solar Imaging in Ca II K and Red Continuum Filters',size=16)
plt.show()

#-----------------------------------------------------------------------
print "End of Problem A(i)"
print "Press enter to begin Problem A(ii)"
my.pause()
#-----------------------------------------------------------------------
#this section cuts through RDC and Contrast images displaying the intensity values along those cuts

mid = np.shape(specTime(1702)[0])[0]/2 #finds the middle of the data

plt.figure(2) #this figure displays both images created
for i in range(2):
	plt.subplot(211+i)
	plt.plot(specTime(1702)[i][mid])
	plt.plot(specTime(1702)[i+2][mid])
	plt.xlim(0,2048)
	plt.xlabel('Pixels')
plt.suptitle('Horizontal Cuts through RDC (top) and Contrast (bottom) Data',size=16)
plt.show()

'''
David: Remember to incude a y label.
'''

#-----------------------------------------------------------------------
print "End of Problem A(ii)"
print "Press enter to begin Problem B(i)"
my.pause()
#-----------------------------------------------------------------------
#this sections takes specific header data from the fits files

def getHeader(time,color,typ): #this function pulls the header from the desired file
	headerData = fits.getheader("20050702."+str(time)+".HW."+color+".P."+typ+".fits")
	return headerData

times = [1702,1720,1740,1802,1820,1840] #this is a list of all 6 times available in given files
headerDataR = []
headerDataK = []
for i in times:
	headerDataR.append(getHeader(i,'R','rdc'))
	headerDataK.append(getHeader(i,'K','rdc'))

Rvals = [[],[],[]]
Kvals = [[],[],[]]
for i in range(6):
	Rvals[0].append(headerDataR[i]['JULDATE']) 
	Rvals[1].append(headerDataR[i]['AVGWIDTH']) 
	Rvals[2].append(headerDataR[i]['RMSA']) 
	Kvals[0].append(headerDataK[i]['JULDATE'])
	Kvals[1].append(headerDataK[i]['AVGWIDTH'])
	Kvals[2].append(headerDataK[i]['RMSA'])

#-----------------------------------------------------------------------
print "End of Problem B(i)"
print "Press enter to begin Problem B(ii)"
my.pause()
#-----------------------------------------------------------------------
#this section takes the data found in the previous section and plots it

plt.figure(3) #this figure contains both plots of compared data
plt.subplot(211)
plt.suptitle('Comparison of Header Values (Red = red line, Ca II = blue line)',size=16)
plt.plot(Rvals[0],Rvals[1],color='r')
plt.plot(Kvals[0],Kvals[1],color='b')
plt.scatter(Rvals[0],Rvals[1],color='r',marker='s',s=30)
plt.scatter(Kvals[0],Kvals[1],color='b',marker='s',s=30)
my.basicPlot('JULDATE vs AVGWIDTH','','AVGWIDTH',False)

plt.subplot(212)
plt.plot(Rvals[0],Rvals[2],color='r')
plt.plot(Kvals[0],Kvals[2],color='b')
plt.scatter(Rvals[0],Rvals[2],color='r',marker='s',s=30)
plt.scatter(Kvals[0],Kvals[2],color='b',marker='s',s=30)
my.basicPlot('JULDATE vs RMSA','JULDATE','RMSA',True)

#-----------------------------------------------------------------------
print "End of Problem B(ii)"
print "Press enter to begin Problem B(iii)"
my.pause()
#-----------------------------------------------------------------------
#this section takes actual center of first three red rdc files and compares sharpness of their edges

redCenter = []
for i in range(3): #this loop pulls the center values of first three red files
	redCenter.append(headerDataR[i]['CENTER_Y'])

plt.figure(4) #this figure plots these intensities at the found center
for i in range(3):
	plotData = specTime(times[i])[0]
	plt.subplot(211)
	plt.plot(plotData[redCenter[i]])
	plt.xlim(1920,1950)
	plt.ylim(0,2000)
	my.basicPlot('','Pixels','Intensity',False)
	
	plt.subplot(212)
	plt.plot(plotData[redCenter[i]])
	plt.xlim(120,150)
	plt.ylim(0,2000)
	my.basicPlot('','Pixels','Intensity',False)
plt.suptitle('Image Quality of Intensity Focused on Edges',size=16)
plt.show()

#-----------------------------------------------------------------------
print "End of Problem B(iii)"
print "Press enter to begin Problem C(i)"
my.pause()
#-----------------------------------------------------------------------
#this section takes the contrast images and plots histograms of variations in intensity accross whole star

def makeHist(data,n): #this funciton gets necessary hist data while ignoring the dark edges of fits data
	hist,binv = np.histogram(data[np.where(data > np.min(data))],bins = n,normed = True)
	barwidth = binv[1]-binv[0]
	return hist,binv,barwidth

def removeDarks(data):
	return data[np.where(data>np.min(data))]

'''
David: You can achieve the same as above using
return data[data>np.min(data)]
The np.where function is useful only if you need to assign its result to a variable.
'''

#specTime(time) returns 2-D arrays, so flatten() is used to turn it into 1-D (easier to work with)
nbins = 100

'''
David: A rule of thumb for making histograms is to make the bins as small as you can while still
capturing the shape of the data. Since you have over 2 million counts, I think you could use
even finer bins.
'''

conR = [specTime(1702)[1].flatten(),specTime(1720)[1].flatten(),specTime(1740)[1].flatten(),specTime(1802)[1].flatten(),specTime(1820)[1].flatten(),specTime(1840)[1].flatten()]
conK = [specTime(1702)[3].flatten(),specTime(1720)[3].flatten(),specTime(1740)[3].flatten(),specTime(1802)[3].flatten(),specTime(1820)[3].flatten(),specTime(1840)[3].flatten()]



plt.figure(5) #this figure displays both of the histograms created of the first three red and ca ii files

plt.subplot(211)
plt.suptitle('Distribution of Intensity from Contrast Image Data',size=16)
plt.hist(removeDarks(conR[0]),bins=nbins,histtype='step',normed=True,range=(-0.04,0.04))
plt.hist(removeDarks(conR[1]),bins=nbins,histtype='step',normed=True,range=(-0.04,0.04))
plt.hist(removeDarks(conR[2]),bins=nbins,histtype='step',normed=True,range=(-0.04,0.04))
plt.xlim(-0.04,0.04)
my.basicPlot('','Intensity','Amount per Bin',False)
plt.subplot(212)
plt.hist(removeDarks(conK[0]),bins=nbins,histtype='step',normed=True,range=(-0.07,0.25))
plt.hist(removeDarks(conK[1]),bins=nbins,histtype='step',normed=True,range=(-0.07,0.25))
plt.hist(removeDarks(conK[2]),bins=nbins,histtype='step',normed=True,range=(-0.07,0.25))
plt.xlim(-0.07,0.25)
my.basicPlot('','Intensity','Amount per Bin',True)

#-----------------------------------------------------------------------
print "End of Problem C(i)"
print "Press enter to begin Problem C(ii)"
my.pause()
#-----------------------------------------------------------------------
#this sectioin plots the second, third fourth, and fifth moment of the intensity distribution after removing the dark area
momR = [[],[],[],[]]
momK = [[],[],[],[]]

#this loop appends the calculated moment using the scipy function
for i in range(6):
	momR[0].append(sci.moment(removeDarks(conR[i]),2))
	momK[0].append(sci.moment(removeDarks(conK[i]),2))
	momR[1].append(sci.moment(removeDarks(conR[i]),3))
	momK[1].append(sci.moment(removeDarks(conK[i]),3))
	momR[2].append(sci.moment(removeDarks(conR[i]),4))
	momK[2].append(sci.moment(removeDarks(conK[i]),4))
	momR[3].append(sci.moment(removeDarks(conR[i]),5))
	momK[3].append(sci.moment(removeDarks(conK[i]),5))

plt.figure(6) #this figure plots the second moment of the data
plt.suptitle('Second Moment of Normalized PDFs of On-Disk Contrast',size=16)
plt.subplot(211)
plt.plot(Rvals[0],momR[0])
my.basicPlot('','','Second Moment',False)
plt.subplot(212)
plt.plot(Kvals[0],momK[0])
my.basicPlot('','JULDATE','Second Moment',True)

def subPlot(x,y,loc): #this function is used to simplfy plotting of the following figure
	sp = plt.subplot(loc)
	plt.plot(x,y)
	plt.ylim(min(y),max(y))
	sp.yaxis.set_major_locator(ticker.MultipleLocator(base=(max(y)-min(y))/4.0)) #these lines set tick frequency
	sp.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))

'''
David: Include markers to show where data points are.
'''

plt.figure(7) #this figure  plots the other three moments as a function of time
plt.suptitle('Higher Moments of Normalized PDFs of On-Disk Contrast',size=16)
subPlot(Rvals[0],momR[1],321)
my.basicPlot('','','Third Moment',False)
subPlot(Rvals[0],momK[1],322)
my.basicPlot('','','',False)
subPlot(Rvals[0],momR[2],323)
my.basicPlot('','','Fourth Moment',False)
subPlot(Rvals[0],momK[2],324)
my.basicPlot('','','',False)
subPlot(Rvals[0],momR[3],325)
my.basicPlot('','JULDATE','Fifth Moment',False)
subPlot(Rvals[0],momK[3],326)
my.basicPlot('','JULDATE','',True)

#-----------------------------------------------------------------------
print "End of Problem C(ii)"
print "Press enter to begin Problem C(iii)"
my.pause()
#-----------------------------------------------------------------------
#this section plots the second moment as a function of AVGWIDTH
plt.figure(8)
plt.suptitle('Second Moment of Contrast Intensity vs. AVGWIDTH',size=16)
plt.subplot(211)
plt.scatter(Rvals[1],momR[0],s=40)
plt.ylim(min(momR[0])*0.99,max(momR[0])*1.01) #this limit scaling is merely so the plot looks better
my.basicPlot('','','Second Moment',False)
plt.subplot(212)
plt.scatter(Rvals[1],momK[0],s=40)
plt.ylim(min(momK[0])*0.99,max(momK[0])*1.01)
my.basicPlot('','AVGWIDTH','SecondMoment',True)

#-----------------------------------------------------------------------
print "End of Problem C(iii)"
print "Press enter to end project"
my.pause()
#-----------------------------------------------------------------------

'''
Grade Breakdown
Completeness: (27/25)
Code ran from start to finish.
Created all plots.
Calculated second moment.
Provided external function file.
Calculated first four central moments. (+1)
Displayed plots of solar limb. (+1)

Logic and Readability: (10/10)
Second moment calculated correctly.
Code is well commented.
Code is formatted in a readable fashion.
Used descriptive variable names.
Plots were well formatted, but subplots could be labeled more clearly.

Efficiency: (5/5)
Used array operations instead of loops where appropriate.
Placed repetitive tasks in functions.
Only used variables that were needed.

Overall: (42/40)
'''




