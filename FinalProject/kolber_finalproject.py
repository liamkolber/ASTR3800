import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci
import scipy.signal as signal
import scipy.optimize as opt
import random

my.introduction('28 Nov, 2016','14 Dec, 2016')
#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin"
my.pause()
#-----------------------------------------------------------------------

def fourier(flux,wavelength): #this function removes the trend from the data
	fftFlux = np.fft.fft(flux)
	freqs = np.fft.fftfreq(flux.size)
	fftFlux[freqs<0.01] = 0 #removes frequencies greater 0.1 (arbitrary through visual anaylsis)
	newFlux = np.real(np.fft.ifft(fftFlux))
	return newFlux

def gaussian(x,a,x0,sigma): #this function is used for curve fit in fitting a gaussian to the peaks
	return a*np.exp(-(x-x0)**2/(2*sigma**2))

def uncertaintyCalculator(sigZ,z): #this function calculates the redshift and its uncertainty 
	mu = np.sum(z/(sigZ**2))/np.sum(1/(sigZ**2))
	sigMu = np.sqrt(1/np.sum(1/(sigZ**2)))
	return mu,sigMu

def multipleLineCheck(firstIndex): #this function is used for analyzing how uncertainty scales with the number of lines
	newRedshiftArray = [redshift[firstIndex]] #creates an array using the redshift of the given index
	newSigmaArray = [sigma[firstIndex]] #same as above
	resultantSigma = []
	resultantShift = []
	randomIndex = random.sample(xrange(len(redshift)),len(redshift)) #generates a sequence of unique random numbers in given range
	resultantSigma.append(uncertaintyCalculator(np.array(newSigmaArray),np.array(newRedshiftArray))[1])
	#sorry for adding more loops/condtions within loops behind a function but I couldn't think of another way
	k = 0
	for i in randomIndex: #this loop appends the random redshift/sigma and calculates uncertainty before appending another
		if i != firstIndex: #ensures that it's not reccalculating with the first index
			newRedshiftArray.append(redshift[i])
			newSigmaArray.append(sigma[i])
			vals = (uncertaintyCalculator(np.array(newSigmaArray),np.array(newRedshiftArray)))
			resultantShift.append(vals[0])
			resultantSigma.append(vals[1])
	if sorted(resultantSigma,reverse=True) != resultantSigma: #checks if the uncertainties are descending
		print "No clear trend\n"
		return False,resultantSigma,resultantShift
	return True,resultantSigma,resultantShift

#-----------------------------------------------------------------------
print "All necessary function are loaded"
print "Press enter to continue"
my.pause()
#-----------------------------------------------------------------------
#inital data analysis

#read in data
data = fits.getdata('spec-0429-51820-0056.fits')
flux = np.array([i[0] for i in data])
wavelength = np.array([10**i[1] for i in data])
restWavelengths = np.genfromtxt('linelist-0429-51820-0056.csv',delimiter=',',skip_header=1).T[1]

#remove the trend from the data using the fourier function
detrended = abs(fourier(flux,wavelength))

#cross correlate the rest wavelengths with the spectrum to find likely shift
shift = (wavelength[1]-wavelength[0])/wavelength[0]
sums = []
for i in range(1000):
	lambdaNew = restWavelengths + i*shift*restWavelengths
	sums.append(np.sum(np.interp(lambdaNew,wavelength,detrended)))

#use the cross correlation to determine redshift
index = sums.index(max(sums))
z = shift*index
print "The estimated redshift via cross correlation: ",z,'\n'

#plot the data before and after the trend is removed
plt.figure(1,figsize=(10, 8))
plt.subplot(211)
plt.plot(wavelength,flux)
plt.scatter(wavelength,flux,s=5)
[plt.axvline(i,color='red') for i in restWavelengths]
plt.xlim(min(wavelength),max(wavelength))
my.basicPlot('Original Spectral Data','','Flux',False)

plt.subplot(212)
[plt.axvline(i,color='red',label='Rest Wavelengths') for i in restWavelengths]
plt.plot(wavelength,detrended)
plt.scatter(wavelength,detrended,s=5)
[plt.axvline(i+shift*index*i,color='gold') for i in restWavelengths]
plt.xlim(min(wavelength),max(wavelength))
my.basicPlot('Detrended Spectral Data','Wavelength','Flux',False)
plt.tight_layout()
plt.show()

#-----------------------------------------------------------------------
print "All data has been adjusted"
print "Press enter to continue"
my.pause()
#-----------------------------------------------------------------------
#guassian fit around wavelength guesses to get a more precise reading

j=2 #this will be used as a counter for the plots
redshift = []
sigma = []
for i in restWavelengths[2:]: #this loop makes a gaussian fit to the peak with a specified range of data
	guessShift = i*(1+z) #the estimated shifted wavelength

	#the range in which the fit will be analyzed (chosen through visual analysis of most prominent peaks)
	lambdaRange = wavelength[(wavelength>(guessShift-10))*(wavelength<(guessShift+10))]
	lambdaInd = np.where((wavelength>(guessShift-10))*(wavelength<(guessShift+10)))
	fluxRange = detrended[lambdaInd]

	try: #this try is here because sometimes the fit can't be made and the function errors
		parameters = [1,guessShift,10/(2*np.sqrt(2*np.log(2)))] #used for curve_fit()
		popt,pcov = opt.curve_fit(gaussian,lambdaRange,fluxRange,parameters)
		gausShift = popt[1]/i-1 #redshift based on gaussian fit
		gausSig = np.sqrt(pcov[1,1])/i #sigma of gaussian fit
		if gausSig < gausShift/100: #ensures uncertainty is less than 1% of redshift
			redshift.append(gausShift)
			sigma.append(gausSig)
			#this plot is made after the fit so that it doesn't create a plot if the fit returns an error
			plt.figure(j,figsize=(10, 8)) 
			#these ranges are used simply to plot the actaul data behind the fit
			largerRange = wavelength[(wavelength>(guessShift-30))*(wavelength<(guessShift+30))]
			largerInd = np.where((wavelength>(guessShift-30))*(wavelength<(guessShift+30)))
			largerFluxRange = detrended[largerInd]
			plt.plot(largerRange,largerFluxRange,'o:',label='Flux Data')
			plt.axvline(guessShift,color='gold',label='Estimated Shift')
			plt.plot(lambdaRange,gaussian(lambdaRange,*popt),'ro:',label='Gaussian Fit')
			plt.xlim(min(largerRange),max(largerRange))
			my.basicPlot('Prominent Peaks Overlayed with Gaussian Fit','Wavelength','Flux',False)
			plt.tight_layout()
			j+=1 #this counter is used for figure making
	except(RuntimeError):
		pass
	plt.legend(loc='upper right')
plt.show()

uncertainty = uncertaintyCalculator(np.array(sigma),np.array(redshift)) #call to calculate redshift/uncertainty
print '\nRedshift = ',uncertainty[0], ' with an estimated error of: +-',uncertainty[1],'\n'

#-----------------------------------------------------------------------
print "A redshift has been determined"
print "Press enter to continue"
my.pause()
#-----------------------------------------------------------------------
#Determine how the uncertainty scales with the number of spectral lines used

uncertaintyScalingArray = []
checks = []
plt.figure(j,figsize=(10, 8)) #this figre displays the change in uncertainties
for i in range(len(redshift)): #this loop calls multipleLineCheck with different starting values
	checks.append(multipleLineCheck(i))
	plt.plot(np.log10(checks[i][1]))
my.basicPlot('Plot of How the Uncertainty Changes as More Spectrum Lines are Added','Number of Lines','log[uncertainties]',True)

#-----------------------------------------------------------------------
print "The change in uncertainty has been evaluated as number of lines increases"
print "Press enter to continue"
my.pause()
#-----------------------------------------------------------------------
#check for systematic differences

for i in range(len(checks)): #this loop analyzes the chi_squared of the relationship between using different lines
	dz = np.array(checks[i][2]-uncertainty[0])
	sigz1 = np.array(checks[i][1])
	sigz2 = np.array(uncertainty[1])
	chi2 = np.sum(dz**2)/np.sum(sigz1**2+sigz2**2)
	p = my.pval(len(dz),chi2) #probability to exceed
	print 'Chi_squraed: ',chi2,'\n','Probability to Exceed: ',p,'\n'

#-----------------------------------------------------------------------
print "Systematic difference have been analyzed"
print "Press enter to end project"
my.pause()
#-----------------------------------------------------------------------
