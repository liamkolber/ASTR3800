import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci

#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin Problem A(i)"
my.pause()
#-----------------------------------------------------------------------

data = np.genfromtxt('SUN_Velocity.txt')
time = [i[0] for i in data]
time = time - min(time)
vel = [i[1] for i in data]

#-----------------------------------------------------------------------
print "End of Problem A(i)"
print "Press enter to begin Problem A(ii)"
my.pause()
#-----------------------------------------------------------------------

plt.figure(1)
plt.scatter(time,vel,s=5)
my.basicPlot('','','',True)

#-----------------------------------------------------------------------
print "End of Problem A(ii)"
print "Press enter to begin Problem A(iii)"
my.pause()
#-----------------------------------------------------------------------

diff = time - np.roll(time,1)
step = max(diff)
newTime = np.arange(0,max(time),step)
newVel = np.interp(newTime,time,vel)

#-----------------------------------------------------------------------
print "End of Problem A(iii)"
print "Press enter to begin Problem A(iv)"
my.pause()
#-----------------------------------------------------------------------

def timePlots(time,newTime,low,high):
	times = [[],[],[],[]]
	for i in range(len(time)):
		if time[i] >= 11.0 and time[i] < 18.0:
			times[0].append(time[i])
			times[1].append(vel[i])
	for i in range(len(newTime)):
		if newTime[i] >= 11.0 and newTime[i] < 18.0:
			times[2].append(newTime[i])
			times[3].append(newVel[i])
	return times

plt.figure(2)
plt.subplot(221)
plt.scatter(time,vel,s=5)
plt.scatter(newTime,newVel,s=5,color='r')
my.basicPlot('','','',False)

plt.subplot(222)
seven = timePlots(time,newTime,11.0,18.0)
plt.scatter(seven[0],seven[1],s=5)
plt.scatter(seven[2],seven[3],s=5,color='r')
my.basicPlot('','','',False)

plt.subplot(223)
hour6 = (1/24.0)*6
six = timePlots(time,newTime,11.0,11.0+hour6)
plt.scatter(six[0],six[1],s=5)
plt.scatter(six[2],six[3],s=5,color='r')
my.basicPlot('','','',False)

plt.subplot(224)
hour = 1/24.0
first = timePlots(time,newTime,11.0,11.0+hour)
plt.scatter(first[0],first[1],s=5)
plt.scatter(first[2],first[3],s=5,color='r')
my.basicPlot('','','',True)

#-----------------------------------------------------------------------
print "End of Problem A(iv)"
print "Press enter to begin Problem A(v)"
my.pause()
#-----------------------------------------------------------------------





