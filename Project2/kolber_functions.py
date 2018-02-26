import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as sci

def introduction(created,edited):
	
	print "Liam Kolber"
	print "Date created: ", created
	print "Date edited: ", edited
	print "\n"

def basicPlot(ttl,x,y):
	plt.title(ttl)
	plt.xlabel(x)
	plt.ylabel(y)
	plt.show() 

def pause():
	
	print '\r'
	stp = raw_input('Pause')
	print '\r'

def timer():
	print '\r'
	t = input('How long of a timer? ')
	print '\r'
	time.sleep(t)

def pval(nbins, chi): #calculate probability to exceed (number of bins, chi_squared)
	p = 1 - sci.chi2.cdf(x = chi, df = nbins - 1)
	return p

