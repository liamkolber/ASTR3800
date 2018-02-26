import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci
import scipy.signal as signal

my.introduction('7 Nov, 2016','28 Nov, 2016')
#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin Section A"
my.pause()
#-----------------------------------------------------------------------
#this section reads in the raw tranist data and plots it, then it
#	separates the individual transits and plots them on top of each
#	other for easy comparison

time,flux,error = np.genfromtxt('wasp_4b.tsv',delimiter=';',dtype=float,skip_header=36).T

plt.figure(1) #this figure plots the raw data with no manipulation
plt.scatter(time,flux)
my.basicPlot('Raw Data','Julian Time','Flux',True)

#these next few lines attempt to find the times of the transits by
#	looking for significant gaps in time.
dt = time - np.roll(time,1)
transitStart = time[dt>=1]
transitStart = np.append(np.min(time),transitStart)
transitEnd = np.roll(time,1)[dt>=1]
transitEnd = np.append(transitEnd,np.max(time))

plt.figure(2) #this figure plots the transits after isolating times 
c = ['blue','green','red','gold','cyan','black']
rawFlux = []
times = []
for i in range(6):
	times.append(time[(time>=transitStart[i])*(time<=transitEnd[i])] - np.min(time[(time>=transitStart[i])])) #the minimum is subtracted to shift times to same start
	rawFlux.append(flux[(time>=transitStart[i])*(time<=transitEnd[i])])
	plt.scatter(times[i],rawFlux[i],color = c[i],label=i+1)
plt.legend(loc='lower right')
my.basicPlot('Transits Adjusted','Julian Time','Flux',True)

#-----------------------------------------------------------------------
print "End of Section A"
print "Press enter to begin Section B"
my.pause()
#-----------------------------------------------------------------------
#this section attempts to allign the transits through cross correlation
#	with two separate methods: physical space and spectral space

#the second transit is used as the reference data, so the interpolation
#	and cross correlation are done relative to that data set
secondTime = time[(time>=transitStart[1])*(time<=transitEnd[1])] - np.min(time[(time>=transitStart[1])])
dt = np.median(secondTime - np.roll(secondTime,1))
interpTime = np.arange(0,max(secondTime),dt) #make a uniform time distribution based on second transit
interpFlux = []
for i in range(6): #this loop interpolates the flux values of all transits
	ogTime = time[(time>=transitStart[i])*(time<=transitEnd[i])] - np.min(time[(time>=transitStart[i])])
	ogFlux = flux[(time>=transitStart[i])*(time<=transitEnd[i])]
	interpFlux.append(np.interp(interpTime,ogTime,ogFlux))

def crossCorr(interp): #this function rolls the data and multiplies it by the unrolled data and sums the array
	allSums = []
	for i in range(len(secondTime)):
		allSums.append(np.sum(interpFlux[1]*np.roll(interp,-i)))
	return allSums

sums = [crossCorr(interpFlux[0]),crossCorr(interpFlux[1]),crossCorr(interpFlux[2]),
	crossCorr(interpFlux[3]),crossCorr(interpFlux[4]),crossCorr(interpFlux[5])]
shifts = dt*np.array([np.where(i==max(i))[0][0] for i in sums]) #this finds the index of the max sum (most optimal shift) and multiplies it by dt in order to find how much (time-wise) to adjust transit data
shiftedTimes = [interpTime-i for i in shifts]

plt.figure(3) #this figure plots the transits again but now with the shifted times so line up
c = ['blue','green','red','gold','cyan','black']
for i in range(6):
	plt.scatter(shiftedTimes[i],interpFlux[i],color=c[i],label=i+1)
plt.legend(loc='lower right')
my.basicPlot('Time Shift Using Physical Space Cross Correlation','Julian Time','Flux',False)

corrThm = []
for i in interpFlux: #this loop attempts the cross correlation using a different method, but time adjustment method is the same as before
	compared = np.fft.fft(i)
	reference = np.fft.fft(interpFlux[1])
	corrThm.append(np.real(np.fft.ifft(compared*np.conjugate(reference))))
CCshifts = dt*np.array([np.where(i==max(i))[0][0] for i in corrThm])
CCshiftedTimes = [interpTime-i for i in CCshifts]

plt.figure(4) #this figure plots the cross correlation data to visually see where the shift that needs to be applied
c = ['blue','green','red','gold','cyan','black']
for i in range(6):
	plt.scatter(interpTime,corrThm[i],color=c[i],label=i+1)
plt.legend(loc='lower right')
my.basicPlot('Cross Correlation','Julian Time','Flux',False)

plt.figure(5) #this figure applies the shift to the data and overplots each transit
c = ['blue','green','red','gold','cyan','black']
for i in range(6):
	plt.scatter(CCshiftedTimes[i],interpFlux[i],color=c[i],label=i+1)
plt.legend(loc='lower right')
my.basicPlot('Time Shift Using Spectral Space Cross Correlation','Julian Time','Flux',True)

#-----------------------------------------------------------------------
print "End of Section B"
print "Press enter to begin Section C"
my.pause()
#-----------------------------------------------------------------------
#this section analyzes the error in the transit data

def fluxSeparater(data,time): #this function separates the non-transit data from the transit data to analyze separately later
	#calculate the left time limit
	leftT = data[(time<=0.01)]
	meanLeft = np.mean(leftT)
	varianceLeft = np.sum((leftT-meanLeft)**2)/len(leftT)
	sigmaLeft = np.sqrt(varianceLeft)
	usableNTLeft = data[(data<=(meanLeft+sigmaLeft))*(data>=(meanLeft-sigmaLeft))*(time<max(time)/2)]
	leftTimeLimit = max(time[(data==usableNTLeft[len(usableNTLeft)-1])*(time<max(time)/2)])

	#calculate the right time limit
	rightT = data[(time>=(max(time)-0.01))]
	meanRight = np.mean(rightT)
	varianceRight = np.sum((rightT-meanRight)**2)/len(rightT)
	sigmaRight = np.sqrt(varianceRight)
	usableNTRight = data[(data<=(meanRight+sigmaRight))*(data>=(meanRight-sigmaRight))*(time>max(time)/2)]
	rightTimeLimit = min(time[(data==usableNTRight[0])*(time>max(time)/2)])

	#calculate the limits for the transit
	transitTime = time[(time>=leftTimeLimit)*(time<=rightTimeLimit)]
	center = np.median(transitTime)
	spread = 0.005
	smallRange = time[(time>=(center-spread))*(time<=center+spread)]
	while len(smallRange)<20:
		spread += 0.001
		smallRange = time[(time>=(center-spread))*(time<=(center+spread))]
	transitData = data[(time>=min(smallRange))*(time<=max(smallRange))]
	meanTransit = np.mean(transitData)
	varianceTransit = np.sum((transitData-meanTransit)**2)/len(transitData)
	sigmaTransit = np.sqrt(varianceTransit)
	usableTransit = data[(data<=(meanTransit+sigmaTransit))*(data>=(meanTransit-sigmaTransit))]
	rightTransitTime = max(time[(data==usableTransit[len(usableTransit)-1])])
	leftTransitTime = min(time[(data==usableTransit[0])])

	#calculate the ratio
	nonTransitFlux = data[(time<=leftTimeLimit)+(time>=rightTimeLimit)]
	transitFlux = data[(time>=leftTransitTime)*(time<=rightTransitTime)]
	ratio = np.mean(transitFlux)/np.mean(nonTransitFlux)

	#return flux data for sigma calculation and means for radius calculation
	return transitFlux,nonTransitFlux,ratio,np.mean(transitFlux),np.mean(nonTransitFlux)

def uncertaintyCalculator(t,non): #this function calculates the uncertainty in the ratios found above (non=nontransit, t=transit)
	#nontransit sigma
	Nnon = len(non)*1.0
	mnon = (1/Nnon)*(np.sum(non))
	variancenon = (np.sum(non**2)*(1/(Nnon-1)) - mnon**2)
	sigMnon = np.sqrt(variancenon)/np.sqrt(Nnon)

	#tansit sigma
	Nt = len(t)*1.0
	mt = (1/Nt)*(np.sum(t))
	variancet = (np.sum(t**2)*(1/(Nt-1)) - mt**2)
	sigMt = np.sqrt(variancet)/np.sqrt(Nt)

	dratio_dnon = np.mean(t)
	dratio_dt = 1/np.mean(non)
	ratio = np.sqrt((dratio_dt*sigMt)**2 + (dratio_dnon*sigMnon)**2)

	return sigMt,sigMnon,ratio

def radiusDensityCalculator(It,sigIt,Inon,sigInon,number):
	SolarRadius = 7*(10**8) #m
	rs = 0.87*SolarRadius
	sigrs = 0.04*SolarRadius
	SolarMass = 2*(10**30) #kg
	JupiterMass = 2*(10**27) #kg
	ms = 0.85*SolarMass
	sigms = 0.11*SolarMass
	mp = 1.21*JupiterMass #mass of WASP-4b
	sigmp = 0.12*JupiterMass

	rp = np.sqrt((rs**2)*(1-It/Inon))
	#partial derivatives for sigma calculation:
	drp_drs = np.sqrt(1-It/Inon)
	drp_dIt = -rs/(2*Inon*drp_drs)
	drp_dInon = (-rs*It)/(2*(Inon**2)*drp_drs)
	#put it all together
	sigrp = np.sqrt((drp_drs*sigrs)**2 + (drp_dIt*sigIt)**2 + (drp_dInon*sigInon)**2)

	#now the density
	dens = mp/((4/3)*(np.pi)*(rp**3))
	#and the partials for uncertainty
	ddens_dmp = 1/((4/3)*(np.pi)*(rp**3))
	ddens_drp = -3/((4/3)*(np.pi)*(rp**4))
	#put it together
	sigdens = np.sqrt((ddens_dmp*sigmp)**2 + (ddens_drp*sigrp)**2)

	if number==5:
		print "Combined Transit Data:"
	else:
		print "Transit ",number+2,":"
	print "Radius of WASP-4b: ",rp," +/- ",sigrp
	print "Density of WASP-4b: ",dens," +/- ",sigdens,'\n'


	return rp, sigrp, dens, sigdens

combinedData = []
for i in range(6):
	combinedData = np.append(combinedData,interpFlux[i])
combinedDataTime = []
for i in range(6):
	combinedDataTime = np.append(combinedDataTime,shiftedTimes[i]) #original flux and original time vlaues

#ignoring the first transit as it is missing vital data
#calculate the ratio of the transit flux to the nontransit flux
ratios = []
[ratios.append(fluxSeparater(rawFlux[i],times[i])) for i in range(1,6)] #original flux and original time vlaues
ratios.append(fluxSeparater(combinedData,combinedDataTime))

#calculate the ratio of the respective uncertainties
sigRatios = []
[sigRatios.append(uncertaintyCalculator(ratios[i][0],ratios[i][1])) for i in range(6)]

#time to calculate the radius
radiusDensityInfo = []
[radiusDensityInfo.append(radiusDensityCalculator(ratios[i][3],sigRatios[i][0],ratios[i][4],sigRatios[i][1],i)) for i in range(6)]

#-----------------------------------------------------------------------
print "End of Section C"
print "Press enter to end the project"
my.pause()
#-----------------------------------------------------------------------