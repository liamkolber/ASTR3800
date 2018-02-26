import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as sci

def introduction(created,edited): #display editing information
	print "Author: Liam Kolber"
	print "Date created: ",created
	print "Date edited: ",edited, '\n'

def basicPlot(ttl,x,y,s): #simple plot labeling function
	plt.title(ttl)
	plt.xlabel(x)
	plt.ylabel(y)
	if s == True:
		plt.show() 

def pause(): #pause the program and ait for user input to continue
	print '\r'
	stp = raw_input('Pause')
	print '\r'

def timer(): #basic timer
	print '\r'
	t = input('How long of a timer? ')
	print '\r'
	time.sleep(t)

def pval(nbins, chi): #calculate probability to exceed (number of bins, chi_squared)
	p = 1 - sci.chi2.cdf(x = chi, df = nbins - 1)
	return p

