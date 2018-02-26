import kolber_functions as my
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats as sci

my.introduction('5 Oct, 2016','TBD')
#-----------------------------------------------------------------------
print "Beginning Project"
print "Press enter to begin Problem A(i) and A(ii) and A(iii)"
my.pause()
#-----------------------------------------------------------------------
#this section reads in the data and produces a plot displaying the daat with a line of best fit

def bestFit(x,y): #this function calculates the values necessary to create a bist fit line
	N = len(x)*1.0
	S1 = N
	Sx = np.sum(np.array(x))
	Sy = np.sum(np.array(y))
	Sxy = np.sum(np.array(x)*np.array(y))
	Sxx = np.sum(np.array(x)**2)

	a = (Sy*Sxx-Sx*Sxy)/(S1*Sxx-Sx**2)
	b = (S1*Sxy-Sx*Sy)/(S1*Sxx-Sx**2)

	eq = b*(np.array(x))+a
	dy = y-eq
	sig = np.sqrt(sci.moment(dy,2)*(N/(N-1)))

	'''
	David: Technically, you should divide by N-2 instead of N-1, but it isn't a big difference.
	'''

	sigS = [S1/(sig**2),Sx/(sig**2),Sy/(sig**2),Sxy/(sig**2),Sxx/(sig**2)]

	sigB = np.sqrt((sigS[0]/(sigS[0]*sigS[4]-sigS[1]**2)))
	sigA = np.sqrt((sigS[4]/(sigS[0]*sigS[4]-sigS[1]**2)))
	
	return a,sigA,b,sigB,eq,sig

def bestFitPlot(x,y,c):#this function takes the datat from bestFit and plots it properly
	a,sigA,b,sigB,eq,sig = bestFit(x,y)
	
	plt.scatter(x,y,label='Original data',color=c)
	plt.plot(x,eq,c,label='Fitted line')
	upper = (b+sigB)*(np.array(x))+(a+sigA)
	lower = (b-sigB)*(np.array(x))+(a-sigA)
	plt.plot(x,upper,c,linestyle='dashed')
	plt.plot(x,lower,c,linestyle='dashed',label='68% Limits')
	return eq,sig

'''
David: When plotting best fit lines and 68% confidence lines, it is more efficient to just use
2 points per line rather than calculating the y-value for every x-value.
'''

csvData = np.genfromtxt('hubbleoriginal.csv', delimiter = ",", skip_header = 1)
hubbleOG = [[],[],[],[]]
for i in csvData:
	hubbleOG[0].append(i[1])
	hubbleOG[1].append(i[2])
	hubbleOG[2].append(i[3])
	hubbleOG[3].append(i[4])

plt.figure(1) #this figure displays hubble's raw data
eq,sig = bestFitPlot(hubbleOG[0],hubbleOG[1],'g')
plt.legend(loc='lower right')
my.basicPlot("Hubble's Original Data",'Distance [Mpc]','Recessional Velocity [km/s]',True)

'''
David: Nice use of functions here.
'''

#-----------------------------------------------------------------------
print "End of Problem A(i) and A(ii) and A(iii)"
print "Press enter to begin Problem A(iv) and A(v)"
my.pause()
#-----------------------------------------------------------------------
#this section calculates a chi-squared value and pearson correlation for the data above

def Xsq(y,mu,sig): #calculate X_squared(expected,observed)
	chi = np.sum(((y-mu)**2)/(sig**2))
	p = 1 - sci.chi2.cdf(x = chi, df = len(y) - 2)
	return chi,p

chi,p = Xsq(hubbleOG[1],eq,sig)

print "Chi-squared ",chi
print "Probability to exceed ",p,"\n"

pearson = sci.pearsonr(hubbleOG[0],hubbleOG[1])
print "Pearson correlation: ", pearson[0], "\n"

#-----------------------------------------------------------------------
print "End of Problem A(iv) and A(v)"
print "Press enter to begin Problem A(vi)"
my.pause()
#-----------------------------------------------------------------------

#extra

#-----------------------------------------------------------------------
print "End of Problem A(Vi)"
print "Press enter to begin Problem B(i) and B(ii) and B(iii)"
my.pause()
#-----------------------------------------------------------------------
#this section performs the same tasks as the first section but with the modern data and over plots for comparison

less8c = np.genfromtxt('ned1dlevel5.csv', delimiter = ",", skip_header = 2,dtype=None)
nameLess = np.array([i[0] for i in less8c])
distLess = np.array([i[3] for i in less8c])
vlcLess = np.array([i[11] for i in less8c])
gals = np.unique(nameLess) 
plotDistL = np.zeros(len(gals))
plotVlcL = np.zeros(len(gals))
for i in range(len(gals)):#this averages data for galaxies so that those don't recieve more weight during calculation
	check = nameLess==gals[i]
	plotDistL[i] = np.mean(distLess[check])
	plotVlcL[i] = np.mean(vlcLess[check])

great8c = np.genfromtxt('ned4dlevel5.csv', delimiter = ",", skip_header = 2,dtype=None)
nameGreat = np.array([i[0] for i in great8c])
distGreat = np.array([i[3] for i in great8c])
vlcGreat = np.array([i[10] for i in great8c])
gals = np.unique(nameGreat)
plotDistG = np.zeros(len(gals))
plotVlcG = np.zeros(len(gals))
simpleGreatData = np.zeros(len(gals))
for i in range(len(gals)):
	check = nameGreat==gals[i]
	plotDistG[i] = np.mean(distGreat[check])
	plotVlcG[i] = np.mean(vlcGreat[check])


plt.figure(2) #this figure overplots these raw data points to compare the three sets of data
bestFitPlot(plotDistL,plotVlcL,'g')
bestFitPlot(plotDistG,plotVlcG,'y')
bestFitPlot(hubbleOG[0],hubbleOG[1],'r')
plt.legend(loc='lower right')
my.basicPlot('Recessional Velocities of Galaxies','Distance [Mpc]','Recessional Velocity [km/s]',True)

'''
David: The best-fit lines for the low-z and Hubble data are not visible in the plot because the x-range
of the lines are the same range as the data. The lines would be visible if you used a larger range.

The white spaces at the edge of the plot should be removed using xlim and ylim.
'''

#-----------------------------------------------------------------------
print "End of Problem B(i) and B(ii) and B(iii)"
print "Press enter to begin Problem B(iv) and B(v) nd B(vi) and B(vii)"
my.pause()
#-----------------------------------------------------------------------
#this section calculates the effective velocities of the data

def dopF(v): #this function uses the doppler formula to calculate the recessional velocities of the galaxies
	c=300000.0 #[km/s]
	v=np.array(v)
	rtrn = np.sqrt(1+v/c)/np.sqrt(1-v/c)-1
	return c*rtrn

veffH = dopF(hubbleOG[1])
veffL = dopF(plotVlcL)
veffG = dopF(plotVlcG)

plt.figure(3) #this overplots the effective velocities of the three data sets
Leq,Lsig = bestFitPlot(plotDistL,veffL,'r')
chiL,pL = Xsq(veffL,Leq,Lsig)
Geq,Gsig = bestFitPlot(plotDistG,veffG,'g')
chiG,pG = Xsq(veffG,Geq,Gsig)
Heq,Hsig = bestFitPlot(hubbleOG[0],veffH,'y')
chiH,pH = Xsq(veffH,Heq,Hsig)
plt.legend(loc='lower right')
my.basicPlot('Recessional Velocities Adjusted for Redshift','Distance [Mpc]','Recessional Velocity [km/s]',True)

pearsonL = sci.pearsonr(plotDistL,veffL)
pearsonG = sci.pearsonr(plotDistG,veffG)
pearsonH = sci.pearsonr(hubbleOG[0],veffH)

print 'Less than 1/8c values: \n','Chi-squared: ',chiL,'\nProbability to exceed: ',pL,'\nPearson Coefficient: ',pearsonL[0],'\n'
print 'Greater 1/8c values: \n','Chi-squared: ',chiG,'\nProbability to exceed: ',pG,'\nPearson Coefficient: ',pearsonG[0],'\n'
print "Hubble's original data values: \n",'Chi-squared: ',chiH,'\nProbability to exceed: ',pH,'\nPearson Coefficient: ',pearsonH[0],'\n'

#-----------------------------------------------------------------------
print "End of Problem B(iv) and B(v) and B(vi) and B(vii)"
print "Press enter to begin Problem B(viii) and B(ix)"
my.pause()
#-----------------------------------------------------------------------

#extra

#-----------------------------------------------------------------------
print "End of Problem B(viii) and B(ix)"
print "Press enter to begin Problem C(i-v)"
my.pause()
#-----------------------------------------------------------------------
#this section attempt to correct the data for solar motion and find hubble's constant

ra = np.array(hubbleOG[2])*(np.pi/12)
dec = np.array(hubbleOG[3])*(np.pi/180)

f = np.array([hubbleOG[0],np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
A = np.zeros((4,4))
b = np.zeros(4).T

#these loops create the necessary matrices for the linear algebra calculations
for i in range(16):
	A[i/4,i%4]=np.sum(f[i/4]*f[i%4])
for i in range(4):
	b[i%4]=np.sum(veffH*f[i%4])

plt.figure(4) #this figure plots Hubble's data after making the correction
[H0,X,Y,Z] = np.dot(np.linalg.inv(A),b)
corrV = veffH - (X*np.sin(ra)*np.cos(dec)+Y*np.sin(ra)*np.cos(dec)+Z*np.sin(dec))
corrEq,corrSig = bestFitPlot(hubbleOG[0],corrV,'b')
plt.legend(loc='lower right')
my.basicPlot('Recessional Velocity Corrected for Solar Motion','Distance [Mpc]','Recessional Velocity (corrected) [km/s]',True)

#the chi-squared of the corrected data to check for strength in relationship
chiCorr,pCorr = Xsq(corrV,corrEq,corrSig)
print 'Corrected velocity values: \n','Chi-squared: ',chiCorr,'\nProbability to exceed: ',pCorr,'\n'

'''
David: Because your plotting function was combined with your linear fit function, you could not make a plot
without redoing the fit. For this part, though, we just wanted you to plot corrV against H0*r (without an
intercept).

In general, I like the way you used functions in this project. However, I think in a couple instances having
those functions limited your ability customize different plots and fits. You can get the functions to be more
flexible and do what you want, it might just take a little extra work.
'''

#-----------------------------------------------------------------------
print "End of Problem C(i-v)"
print "Press enter to end Project 4"
my.pause()
#-----------------------------------------------------------------------

'''
Grade Breakdown

Completeness: (25/25)
Code ran from start to finish.
Code created all necessary plots.
Calculated chi-squared and p-values.
Calculated best-fit values for all data sets.

Logic and Readability: (9/10)
Code included descriptive comments.
Variable names were descriptive.
Chi-squared tests performed correctly with correct degrees of freedom.
White spaces need to be removed from edge of plots and best-fit lines for low-z data needs to be extended to higher distances. (-1)

Efficiency: (5/5)
Very good use of functions.
Used array operations instead of loops where appropriate.

Overall (39/40)
'''





