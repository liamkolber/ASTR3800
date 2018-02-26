#import time
import numpy as np
import time

def introduction(created, edited):
	
	print "Liam Kolber"
	print "Date created: ", created
	print "Date edited: ", edited
	print "\n"
	
def pause():
	
	print '\r'
	stp = raw_input('Pause')
	print '\r'

def timer():
	print '\r'
	t = input('How long of a timer? ')
	print '\r'
	time.sleep(t)

