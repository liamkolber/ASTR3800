import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci
import scipy.interpolate as inter

#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin Problem A(i)"
my.pause()
#-----------------------------------------------------------------------
#load in data and shift time to start at zero
time,vel = np.genfromtxt('SUN_Velocity.txt').T
time = time - min(time)

#-----------------------------------------------------------------------
print "End of Problem A(i)"
print "Press enter to begin Problem A(ii)"
my.pause()
#-----------------------------------------------------------------------
#plot the data attained from the previous section
plt.figure(1)
plt.plot(time,vel)
my.basicPlot('Doppler Velocity of Solar Surface','Time [s]','Velocity [m/s]',True)

#-----------------------------------------------------------------------
print "End of Problem A(ii)"
print "Press enter to begin Problem A(iii)"
my.pause()
#-----------------------------------------------------------------------
#establish a time-step by finding the distance between each point and finding the max
diff = time - np.roll(time,1)
step = max(diff)
newTime = np.arange(0,max(time),step) #new time array with new time-step
newVel = np.interp(newTime,time,vel)

#-----------------------------------------------------------------------
print "End of Problem A(iii)"
print "Press enter to begin Problem A(iv)"
my.pause()
#-----------------------------------------------------------------------
#plot various time ranges to see the spectrum at different scales

def AivPlots(VEL,upperLim): #this function is just to simplify plot making
	plt.plot(time[(time>=11.0)*(time<=upperLim)],vel[(time>=11.0)*(time<=upperLim)])
	plt.plot(newTime[(newTime>=11.0)*(newTime<=upperLim)],VEL[(newTime>=11.0)*(newTime<=upperLim)],color='r')

plt.figure(2)
plt.suptitle('Doppler Velocity of Solar Surface with Adjusted Time Step',fontsize=20)
plt.subplot(411)
plt.plot(time,vel)
plt.plot(newTime,newVel,color='r')
my.basicPlot('Full Time Series','','Velocity [m/s]',False)

plt.subplot(412)
AivPlots(newVel,18.0)
my.basicPlot('Days 11-18','','Velocity [m/s]',False)

plt.subplot(413)
hour6 = (1/24.0)*6
AivPlots(newVel,11.0+hour6)
my.basicPlot('Day 11 (First 6 Hours)','','Velocity [m/s]',False)

plt.subplot(414)
hour = 1/24.0
AivPlots(newVel,11.0+hour)
my.basicPlot('Day 11 (First Hour)','Time [s]','Velocity [m/s]',True)

#-----------------------------------------------------------------------
print "End of Problem A(iv)"
print "Press enter to begin Problem A(v)"
my.pause()
#-----------------------------------------------------------------------
#plot various time ranges but using the nearest neighbor method of interpolation 
nearVel = inter.griddata(time,vel,newTime,'nearest')

plt.figure(3)
plt.suptitle('Doppler Velocity of Solar Surface with Nearest Neightbor Sampling',fontsize=20)
plt.subplot(411)
plt.plot(time,vel)
plt.plot(newTime,newVel,color='r')
my.basicPlot('Full Time Series','','Velocity [m/s]',False)

plt.subplot(412)
AivPlots(nearVel,18.0)
my.basicPlot('Days 11-18','','Velocity [m/s]',False)

plt.subplot(413)
hour6 = (1/24.0)*6
AivPlots(nearVel,11.0+hour6)
my.basicPlot('Day 11 (First 6 Hours)','','Velocity [m/s]',False)

plt.subplot(414)
hour = 1/24.0
AivPlots(nearVel,11.0+hour)
my.basicPlot('Day 11 (First Hour)','Time [s]','Velocity [m/s]',True)

#-----------------------------------------------------------------------
print "End of Problem A(v)"
print "Press enter to begin Problem A(vi)"
my.pause()
#-----------------------------------------------------------------------
#remove any trends in the data so that it is centered around zero

fullSeriesFitVals = my.bestFit(newTime,newVel)
fullSeries = newVel-fullSeriesFitVals[4]

sevSeriesFitVals = my.bestFit(newTime[(newTime>=11.0)*(newTime<=18.0)],newVel[(newTime>=11.0)*(newTime<=18.0)])
sevenSeries = newVel[(newTime>=11.0)*(newTime<=18.0)]-sevSeriesFitVals[4]

plt.figure(4)
plt.suptitle('Time/Velocity Data Series with Trend Removed')
plt.subplot(211)
plt.plot(newTime,fullSeries)
my.basicPlot('Full Time Series','','Velocity [m/s]',False)
plt.subplot(212)
plt.plot(newTime[(newTime>=11.0)*(newTime<=18.0)],sevenSeries)
my.basicPlot('Days 11 - 18','Time [s]','Velocity [m/s]',True)

#-----------------------------------------------------------------------
print "End of Problem A(vi)"
print "Press enter to begin Problem B(i)"
my.pause()
#-----------------------------------------------------------------------
#calculate the power spectra of the data sets for full time series and seven day time series

def wave(signal,time,dt): #this function calculates the power spectra
	amp = np.fft.fft(signal)/(len(time))
	powr = np.real(amp*np.conjugate(amp))
	freq = np.fft.fftfreq(powr.size,dt)
	df = (freq[1]-freq[0])
	return freq,powr,df

fullfreq,fullpowr,dfull = wave(fullSeries,newTime,step*(24*3600))
sevfreq,sevpowr,dsev = wave(sevenSeries,newTime[(newTime>=11.0)*(newTime<=18.0)],step*(24*3600))

fullpowr = (fullpowr/np.sum(fullpowr[fullfreq>=0]*dfull))*(fullSeriesFitVals[5]**2)
sevpowr = (sevpowr/np.sum(sevpowr[sevfreq>=0]*dsev))*(sevSeriesFitVals[5]**2)

#-----------------------------------------------------------------------
print "End of Problem B(i)"
print "Press enter to begin Problem B(ii)"
my.pause()
#-----------------------------------------------------------------------
#plot the power spectra in a log-log plot

fullfreq =fullfreq[fullfreq>=0]*1000
sevfreq = sevfreq[sevfreq>=0]*1000

plt.figure(5)
plt.suptitle('Power Spectra of Velocity Time-Series')
plt.subplot(211)
plt.plot(np.log10(fullfreq),np.log10(fullpowr[fullfreq>=0]))
plt.xlim(-4.0,1.0)
plt.ylim(-3.0,6.0)
my.basicPlot('Full Time Series','','log(Power)',False)
plt.subplot(212)
plt.plot(np.log10(sevfreq),np.log10(sevpowr[sevfreq>=0]))
plt.xlim(-4.0,1.0)
plt.ylim(-3.0,6.0)
my.basicPlot('Days 11 - 18','Frequency [mHz]','log(Power)',True)

#-----------------------------------------------------------------------
print "End of Problem B(ii)"
print "Press enter to begin Problem B(iii)"
my.pause()
#-----------------------------------------------------------------------
#plot the p-mode portion of the spectra

plt.figure(6)
plt.suptitle('P-Mode Portion of Spectra')
plt.subplot(211)
plt.plot(np.log10(fullfreq),np.log10(fullpowr[fullfreq>=0]))
plt.xlim(np.log10(2.0),np.log10(7.0))
plt.ylim(-3.0,6.0)
my.basicPlot('Full Time Series','','log(Power)',False)
plt.subplot(212)
plt.plot(np.log10(sevfreq),np.log10(sevpowr[sevfreq>=0]))
plt.xlim(np.log10(2.0),np.log10(7.0))
plt.ylim(-3.0,6.0)
my.basicPlot('Days 11 - 18','Frequency [mHz]','log(Power)',True)

#-----------------------------------------------------------------------
print "End of Problem B(iii)"
print "Press enter to begin Problem B(iv)"
my.pause()
#-----------------------------------------------------------------------
#isolate a small portion of the spectra to focus on the peak

spectraFullFreq = fullfreq[(fullfreq>=3.2)*(fullfreq<=3.4)]
spectraFullPowr = fullpowr[(fullfreq>=3.2)*(fullfreq<=3.4)]
spectraSevFreq = sevfreq[(sevfreq>=3.2)*(sevfreq<=3.4)]
spectraSevPowr = sevpowr[(sevfreq>=3.2)*(sevfreq<=3.4)]

subrangeFullFreq = fullfreq[(fullfreq>=3.25)*(fullfreq<=3.35)]
subrangeFullPowr = fullpowr[(fullfreq>=3.25)*(fullfreq<=3.35)]
subrangeSevFreq = sevfreq[(sevfreq>=3.25)*(sevfreq<=3.35)]
subrangeSevPowr = sevpowr[(sevfreq>=3.25)*(sevfreq<=3.35)]

plt.figure(7)
plt.suptitle('Linear Frequency vs log(Power)')
plt.subplot(211)
plt.plot(spectraFullFreq,np.log10(spectraFullPowr))
plt.plot(subrangeFullFreq,np.log10(subrangeFullPowr),label='Focus Area')
plt.xlim(3.2,3.4)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Full Time Series','','log(Power)',False)
plt.subplot(212)
plt.plot(spectraSevFreq,np.log10(spectraSevPowr))
plt.plot(subrangeSevFreq,np.log10(subrangeSevPowr),label='Focus Area')
plt.xlim(3.2,3.4)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Days 11 - 18','Frequency[mHz]','log(Power)',True)

#-----------------------------------------------------------------------
print "End of Problem B(iv)"
print "Press enter to begin Problem B(vi:A)"
my.pause()
#-----------------------------------------------------------------------
#further isolate a portion to bring even more attention to the peak and immediate surroundings

focusFullFreq = fullfreq[(fullfreq>=3.285)*(fullfreq<=3.315)]
focusFullPowr = fullpowr[(fullfreq>=3.285)*(fullfreq<=3.315)]
focusSevFreq = sevfreq[(sevfreq>=3.285)*(sevfreq<=3.315)]
focusSevPowr = sevpowr[(sevfreq>=3.285)*(sevfreq<=3.315)]
backgroundFullFreq = fullfreq[(fullfreq>=3.25)*(fullfreq<=3.285)+(fullfreq>=3.315)*(fullfreq<=3.35)]
backgroundFullPowr = fullpowr[(fullfreq>=3.25)*(fullfreq<=3.285)+(fullfreq>=3.315)*(fullfreq<=3.35)]
backgroundSevFreq = sevfreq[(sevfreq>=3.25)*(sevfreq<=3.285)+(sevfreq>=3.315)*(sevfreq<=3.35)]
backgroundSevPowr = sevpowr[(sevfreq>=3.25)*(sevfreq<=3.285)+(sevfreq>=3.315)*(sevfreq<=3.35)]

meanBGFull = my.bestFit(backgroundFullFreq,backgroundFullPowr) #array[4] & array[5] are mean and variance respectively
meanBGSev = my.bestFit(backgroundSevFreq,backgroundSevPowr)

plt.figure(8)
plt.suptitle('Focused Range with Highlighted Spectral Peak')
plt.subplot(211)
plt.plot(subrangeFullFreq,np.log10(subrangeFullPowr))
plt.plot(focusFullFreq,np.log10(focusFullPowr),label='Spectral Peak')
plt.xlim(3.25,3.35)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Full Time Series','','log(Power)',False)
plt.subplot(212)
plt.plot(subrangeSevFreq,np.log10(subrangeSevPowr))
plt.plot(focusSevFreq,np.log10(focusSevPowr),label='Spectral Peak')
plt.xlim(3.25,3.35)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Days 11 - 18','','log(Power)',True)
#-----------------------------------------------------------------------
print "End of Problem B(vi:A)"
print "Press enter to begin Problem B(vi:B,C,D)"
my.pause()
#-----------------------------------------------------------------------
#calculate a chi_squared of the peak as compared to its background

chiFull = (np.sum((focusFullPowr-np.mean(backgroundFullPowr))**2))/(meanBGFull[5]**2)
pFull = my.pval(len(focusFullPowr),chiFull)
chiSev = (np.sum((focusSevPowr-np.mean(backgroundSevPowr))**2))/(meanBGSev[5]**2)
pSev = my.pval(len(focusSevPowr),chiSev)

print "Full Time Series:","\n","Chi: ",chiFull,"\n","Probability to exceed: ",pFull,"\n"
print "Seven Day Time Series:","\n","Chi: ",chiSev,"\n","Probability to exceed: ",pSev,"\n"

#-----------------------------------------------------------------------
print "End of Problem B(vi)"
print "Press enter to begin Problem B(vii)"
my.pause()
#-----------------------------------------------------------------------
#repeat the past couple sections with a different focus range to check for other notable features

newFocusFullFreq = fullfreq[(fullfreq>=3.3124)*(fullfreq<=3.325)]
newFocusFullSev = sevfreq[(sevfreq>=3.3124)*(sevfreq<=3.325)]
newFocusFullPowr = fullpowr[(fullfreq>=3.3124)*(fullfreq<=3.325)]
newFocusSevPowr = sevpowr[(sevfreq>=3.3124)*(sevfreq<=3.325)]

plt.figure(9)
plt.suptitle('Focused Range with Highlighted Spectral Peak')
plt.subplot(211)
plt.plot(subrangeFullFreq,np.log10(subrangeFullPowr))
plt.plot(newFocusFullFreq,np.log10(newFocusFullPowr),label='Concerned Area')
plt.xlim(3.25,3.35)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Full Time Series','','log(Power)',False)
plt.subplot(212)
plt.plot(subrangeSevFreq,np.log10(subrangeSevPowr))
plt.plot(newFocusFullSev,np.log10(newFocusSevPowr),label='Concerned Area')
plt.xlim(3.25,3.35)
plt.ylim(-3.0,6.0)
plt.legend(loc='lower right')
my.basicPlot('Days 11 - 18','','log(Power)',False)

newChiFull = (np.sum((newFocusFullPowr-np.mean(newFocusFullPowr))**2))/(meanBGFull[5]**2)
newChiSev = (np.sum((newFocusSevPowr-np.mean(newFocusSevPowr))**2))/(meanBGSev[5]**2)
newpFull = my.pval(len(newFocusFullPowr),newChiFull)
newpSev = my.pval(len(newFocusSevPowr),newChiSev)

print "Full Time Series:","\n","Chi: ",newChiFull,"\n","Probability to exceed: ",newpFull,"\n"
print "Seven Day Time Series:","\n","Chi: ",newChiSev,"\n","Probability to exceed: ",newpSev,"\n"



