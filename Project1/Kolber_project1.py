import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sci
import Kolber_functions as my

print "\nBeginning Project 1\n"
print "Press enter to begin Problem A(i)"
my.pause()

# Problem A(i)
rand1 = []
rand2 = []
rand3 = []
randlist = [rand1, rand2, rand3] #empty array for random numbers
seed1 = 414929
seed2 = 581714
seed3 = 698137
seed = [seed1, seed2, seed3] #array of all three seeds together

#for loop runs through each seed seperately
for j in range (0, 3):
	np.random.seed(seed[j]) #establish seed based on my 3 choices
	limit = 1000 #number of random numbers
	i = 0
	
	#while loop creates the random numbers
	while i < limit :
		num = np.random.uniform(0,1.0)
		randlist[j].append(num)
		print num
		i += 1
	
	print '\n' #to easily differentiate between each seed's list

print "End of Problem A(i)"
print "Press enter to begin Problem A(ii)"
my.pause()

#Problem A(ii)

#using the same seeds as above
jv1 = [] #empty array for each seed
jv2 = []
jv3 = []
JvN_list = [jv1, jv2, jv3] #array of arrays (JvN for John von Neumann)
j = 0

for j in range (0, 3):
	i = 0
	randnum = (seed[j]*1.0)**2 #square the first seed number (multiply by 1.0 to convert to float
	numcheck = seed[j]*1.0 #numcheck saves the seed in case of flatline
	alert = False #this will be used in the check for a flatline
	
	#while loop creates random numbers through JvN method and checks for duplicates/flatline
	while i < 1000:
		mystr = str(randnum**2) #converts square into string to pull middle 6 digits out
		
		#this checks to make sure I can pull 6 digits out (pulling 3:9)
		#if it fails, establishes new seed number
		if len(mystr) < 10:
			newnum = numcheck + 1 #changes the number to be the original seed + 1
			sqnum = (newnum)**2 #runs through math again by squaring
			numcheck += 1 #changes the saved seed to reflect the current, changed seed
		#if the number meets length requirement
		else:
			newnum = int(mystr[3:9]) #pull middle 6 and revert back to int
			sqnum = (newnum)**2 #square those new 6 numbers
		
		#ensures that there are no zeros in the list
		if newnum/1000000.0 == 0:
			newnum = numcheck + 1
			numcheck += 1
		
		JvN_list[j].append(newnum/1000000.0) #append the array to contain the new 6 digit nunbers
		print newnum/1000000.0 #prints the numbers with normalization
		
		#runs through current array to ensure uniqueness
		for k in range (0, len(JvN_list)-1):
			if (k == newnum):
				alert = True
				break
		if alert == True:
			newnum = numcheck + 1 #changes number in case fails test
			numcheck += 1
		
		#prepare for next iteration
		alert = False
		i += 1
		randnum = newnum
		
	print '\n'

print "End of Problem A(ii)"
print "Press enter to begin Problem B(i)"
my.pause()

#the following plots a line graph displaying each number generated
#using numpy.random.uniform
plt.figure(1)
plt.plot(rand1)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using numpy.random.uniform")

plt.figure(2)
plt.plot(rand2)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using numpy.random.uniform")

plt.figure(3)
plt.plot(rand3)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using numpy.random.uniform")

#the following plots a line graph displaying each number generated
#using John von Neumann method
plt.figure(4)
plt.plot(jv1)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using John von Neumann Method")

plt.figure(5)
plt.plot(jv2)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using John von Neumann Method")

plt.figure(6)
plt.plot(jv3)
plt.xlabel("Iteration")
plt.ylabel("Range")
plt.title("Random Number Generation Using John von Neumann Method")

plt.show()

print "End of Problem B(i)"
print "Press enter to begin Problem B(ii)"
my.pause()

#the following plots a histogram to see distribution of generated numbers
#using numpy.random.uniform

plt.figure(7)
plota10 = plt.hist(rand1)
plota40 = plt.hist(rand1, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using numpy.random.uniform")

plt.figure(8)
plotb = plt.hist(rand2)
plotb40 = plt.hist(rand2, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using numpy.random.uniform")

plt.figure(9)
plotc = plt.hist(rand3)
plotc40 = plt.hist(rand3, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using numpy.random.uniform")

#the following plots a histogram to see distribution of generated numbers
#using John von Neumann method
plt.figure(10)
plotd = plt.hist(jv1)
plotd40 = plt.hist(jv1, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using John von Neumann Method")

plt.figure(11)
plote = plt.hist(jv2)
plote40 = plt.hist(jv2, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using John von Neumann Method")

plt.figure(12)
plotf = plt.hist(jv3)
plotf40 = plt.hist(jv3, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("Random Number Generation Using John von Neumann Method")

#the following are combined plots of all 3 runs for each technique
NPRandCom = []
JvNCom = []
NPRandCom.extend(rand1)
NPRandCom.extend(rand2)
NPRandCom.extend(rand3)
JvNCom.extend(jv1)
JvNCom.extend(jv2)
JvNCom.extend(jv3)

plt.figure(13)
plotg = plt.hist(NPRandCom)
plotg40 = plt.hist(NPRandCom, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("RNG Using numpy.random.uniform (all runs combined)")

plt.figure(14)
plotf = plt.hist(JvNCom)
plotf40 = plt.hist(JvNCom, 40, facecolor = 'green')
plt.xlabel("Random Number")
plt.ylabel("Number of Occurances")
plt.title("RNG Using John von Neumann Method (all runs combined)")

plt.show()

print "End of Problem B(ii)"
print "Press enter to begin Problem B(iii,a,i)"
my.pause()

#all of this just to figure out the amount in each bin using numpy.histogram
#from what I saw it returned as a string so I convert it to ints
def removeSpaces(temp): #this function removes the empty spaces in list
	loop = 0
	cnt = 0
	for n in temp:
		if n == '':
			cnt += 1
	while loop < cnt:
		temp.remove('')
		loop += 1

def histSize(n, array): #this function figures out the size of each bin
	bin1 = np.histogram(rand1, bins = n)
	bin2 = np.histogram(rand2, bins = n)
	bin3 = np.histogram(rand3, bins = n)
	bin4 = np.histogram(jv1, bins = n)
	bin5 = np.histogram(jv2, bins = n)
	bin6 = np.histogram(jv3, bins = n)
	bin7 = np.histogram(NPRandCom, bins = n)
	bin8 = np.histogram(JvNCom, bins = n)

	binnum1 = str(bin1[0])
	binnum2 = str(bin2[0])
	binnum3 = str(bin3[0])
	binnum4 = str(bin4[0])
	binnum5 = str(bin5[0])
	binnum6 = str(bin6[0])
	binnum7 = str(bin7[0])
	binnum8 = str(bin8[0])
	binlist = [binnum1,binnum2,binnum3,binnum4,binnum5,binnum6, binnum7, binnum8]
	newlist = []

	for l in binlist:
		
		l = l[1:len(l)-1]
		newlist.append(l)

	for m in newlist:
		temp = m.split(' ')
		removeSpaces(temp)
		array.append(map(int, temp))
		
intlist10 = []
intlist40 = []
intlistHunnit = []
histSize(10, intlist10)
histSize(40, intlist40)
histSize(100,intlistHunnit)

def Xsq(ex, array): #calculate X_squared
	tot = 0
	temp = []
	for p in array:
		tot = ((p - (ex*1.0))**2)/ex
		temp.append(tot)
	Sum = np.sum(temp)
		
	return Sum

def pval(n, chi): #calculate probability to exceed
	p = 1 - sci.chi2.cdf(x = chi, df = n - 1)
	return p

#calculate X-squared for numpy.random
#time to calculate X-squared for 10
[npXsq10_1, npXsq10_2, npXsq10_3] = [Xsq(100, intlist10[0]), Xsq(100, intlist10[1]), Xsq(100, intlist10[2])]
[npXsq40_1, npXsq40_2, npXsq40_3] = [Xsq(25, intlist40[0]), Xsq(25, intlist40[1]), Xsq(25, intlist40[2])]
[npXsq100_1, npXsq100_2, npXsq100_3] = [Xsq(10, intlistHunnit[0]), Xsq(10, intlistHunnit[1]), Xsq(10, intlistHunnit[2])]

#calculate X-squared for JvN
#time to calculate X-squared for 10
[jvXsq10_1, jvXsq10_2, jvXsq10_3] = [Xsq(100, intlist10[3]), Xsq(100, intlist10[4]), Xsq(100, intlist10[5])]
[jvXsq40_1, jvXsq40_2, jvXsq40_3] = [Xsq(25, intlist40[3]), Xsq(25, intlist40[4]), Xsq(25, intlist40[5])]
[jvXsq100_1, jvXsq100_2, jvXsq100_3] = [Xsq(10, intlistHunnit[3]), Xsq(10, intlistHunnit[4]), Xsq(10, intlistHunnit[5])]

[nall10, nall40, nall100] = [Xsq(300, intlist10[6]), Xsq(75, intlist40[6]), Xsq(30, intlistHunnit[6])]
[jall10, jall40, jall100] = [Xsq(300, intlist10[7]), Xsq(75, intlist40[7]), Xsq(30, intlistHunnit[7])]

Xsq_all = [npXsq10_1, npXsq10_2, npXsq10_3, npXsq40_1, npXsq40_2, npXsq40_3, npXsq100_1, npXsq100_2, npXsq100_3, jvXsq10_1, jvXsq10_2, jvXsq10_3, jvXsq40_1, jvXsq40_2, jvXsq40_3, jvXsq100_1, jvXsq100_2, jvXsq100_3, nall10, nall40, nall100, jall10, jall40, jall100]

#calculate P for each bin size of each generated sequence
n10p1 = pval(10, Xsq_all[0])
n10p2 = pval(10, Xsq_all[1])
n10p3 = pval(10, Xsq_all[2])
n40p1 = pval(40, Xsq_all[3])
n40p2 = pval(40, Xsq_all[4])
n40p3 = pval(40, Xsq_all[5])
n100p1 = pval(100, Xsq_all[6])
n100p2 = pval(100, Xsq_all[7])
n100p3 = pval(100, Xsq_all[8])
j10p1 = pval(10, Xsq_all[9])
j10p2 = pval(10, Xsq_all[10])
j10p3 = pval(10, Xsq_all[11])
j40p1 = pval(40, Xsq_all[12])
j40p2 = pval(40, Xsq_all[13])
j40p3 = pval(40, Xsq_all[14])
j100p1 = pval(100, Xsq_all[15])
j100p2 = pval(100, Xsq_all[16])
j100p3 = pval(100, Xsq_all[17])
#calculate P for the combined lists of all three generations
nall10p = pval(10, Xsq_all[18])
nall40p = pval(40, Xsq_all[19])
nall100p = pval(100, Xsq_all[20])
jall10p = pval(10, Xsq_all[21])
jall40p = pval(40, Xsq_all[22])
jall100p = pval(100, Xsq_all[23])

p_all = [n10p1, n10p2, n10p3, n40p1, n40p2, n40p3, n100p1, n100p2, n100p3, j10p1, j10p2, j10p3, j40p1, j40p2, j40p3, j100p1, j100p2, j100p3, nall10p, nall40p, nall100p, jall10p, jall40p, jall100p]

print "End of Problem B(iii,a,i)"
print "Press enter to begin Problem B(iii,a,ii)"
my.pause()

print 'X_squared for each generated sequence:\n', Xsq_all, '\n'
print 'Probability to exceed for each X_squared:\n', p_all, '\n'

print "End of Problem B(iii,a,ii)"
print "Press enter to begin Problem B(iii,b,i)"
my.pause()

#perform serial test by placing generations in bins and checking for
#consecutive and adjacent bin placement (sorry it can look a little confusing)
def serialTest(array):
	b1 = []
	b2 = []
	b3 = []
	b4 = []
	b5 = []
	b6 = []
	b7 = []
	b8 = []
	b9 = []
	b10 = []
	bins = [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10]
	index = 0
	cnt1 = 0
	cnt2 = 0
	cnt3 = 0
	cnt4 = 0
	cnt5 = 0
	cnt6 = 0
	cnt7 = 0
	cnt8 = 0
	cnt9 = 0
	cnt10 = 0
	ncnt1 = 0
	ncnt2 = 0
	ncnt3 = 0
	ncnt4 = 0
	ncnt5 = 0
	ncnt6 = 0
	ncnt7 = 0
	ncnt8 = 0
	ncnt9 = 0
	ncnt10 = 0
	skip = False

	#this loop places numbers in bins and then checks for patterns
	for i in array:
		if i >= 0.0 and i <= 0.1:
			b1.append(i)
			if array[index-1] == b1[len(b1)-2]:
				cnt1 += 1
			if (len(b2)>0) and (array[index-1] == b2[len(b2)-1]):
				ncnt1 += 1
		elif i > 0.1 and i <= 0.2:
			b2.append(i)
			if array[index-1] == b2[len(b2)-2]:
				cnt2 += 1
			if len(b3)>0 and (array[index-1] == b3[len(b3)-1]) and skip == False:
				ncnt2 += 1
				skip == True
			elif len(b1)>0 and (array[index-1] == b1[len(b1)-1]) and skip == False:
				ncnt2 += 1
				skip == True
		elif i > 0.2 and i <= 0.3:
			b3.append(i)
			if array[index-1] == b3[len(b3)-2]:
				cnt3 += 1
			if len(b4)>0 and (array[index-1] == b4[len(b4)-1]) and skip == False:
				ncnt3 += 1
				skip == True
			elif len(b2)>0 and (array[index-1] == b2[len(b2)-1]) and skip == False:
				ncnt3 += 1
				skip == True
		elif i > 0.3 and i <= 0.4:
			b4.append(i)
			if array[index-1] == b4[len(b4)-2]:
				cnt4 += 1
			if len(b5)>0 and (array[index-1] == b5[len(b5)-1]) and skip == False:
				ncnt4 += 1
				skip == True
			elif len(b3)>0 and (array[index-1] == b3[len(b3)-1]) and skip == False:
				ncnt4 += 1
				skip == True
		elif i > 0.4 and i <= 0.5:
			b5.append(i)
			if array[index-1] == b5[len(b5)-2]:
				cnt5 += 1
			if len(b6)>0 and (array[index-1] == b6[len(b6)-1]) and skip == False:
				ncnt5 += 1
				skip == True
			elif len(b4)>0 and (array[index-1] == b4[len(b4)-1]) and skip == False:
				ncnt5 += 1
				skip == True
		elif i > 0.5 and i <= 0.6:
			b6.append(i)
			if array[index-1] == b6[len(b6)-2]:
				cnt6 += 1
			if len(b7)>0 and (array[index-1] == b7[len(b7)-1]) and skip == False:
				ncnt6 += 1
				skip == True
			elif len(b5)>0 and (array[index-1] == b5[len(b5)-1]) and skip == False:
				ncnt6 += 1
				skip == True
		elif i > 0.6 and i <= 0.7:
			b7.append(i)
			if array[index-1] == b7[len(b7)-2]:
				cnt7 += 1
			if len(b8)>0 and (array[index-1] == b8[len(b8)-1]) and skip == False:
				ncnt7 += 1
				skip == True
			elif len(b6)>0 and (array[index-1] == b6[len(b6)-1]) and skip == False:
				ncnt7 += 1
				skip == True
		elif i > 0.7 and i <= 0.8:
			b8.append(i)
			if array[index-1] == b8[len(b8)-2]:
				cnt8 += 1
			if len(b9)>0 and (array[index-1] == b9[len(b9)-1]) and skip == False:
				ncnt8 += 1
				skip == True
			elif len(b7)>0 and (array[index-1] == b7[len(b7)-1]) and skip == False:
				ncnt8 += 1
				skip == True
		elif i > 0.8 and i <= 0.9:
			b9.append(i)
			if array[index-1] == b9[len(b9)-2]:
				cnt9 += 1
			if len(b10)>0 and (array[index-1] == b10[len(b10)-1]) and skip == False:
				ncnt9 += 1
				skip == True
			elif len(b8)>0 and (array[index-1] == b8[len(b8)-1]) and skip == False:
				ncnt9 += 1
				skip == True
		elif i > 0.9:
			b10.append(i)
			if array[index-1] == b10[len(b10)-2]:
				cnt10 += 1
			elif (len(b9)>0) and array[index-1] == b9[len(b9)-1]:
				ncnt10 += 1
				
		index += 1
		skip = False

	Same = [cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,cnt7,cnt8,cnt9,cnt10]
	Next = [ncnt1,ncnt2,ncnt3,ncnt4,ncnt5,ncnt6,ncnt7,ncnt8,ncnt9,ncnt10]
	
	rtrn = [Same, Next]
	return rtrn
	
#call serialTest to find necessary information
[r1Same,r1Next] = serialTest(rand1)
[r2Same,r2Next] = serialTest(rand2)
[r3Same,r3Next] = serialTest(rand3)
[j1Same,j1Next] = serialTest(jv1)
[j2Same,j2Next] = serialTest(jv2)
[j3Same,j3Next] = serialTest(jv3)

#compile results into single lists
same = [r1Same, r2Same, r3Same, j1Same, j2Same, j3Same]
adja = [r1Next, r2Next, r3Next, j1Next, j2Next, j3Next]

print "End of Problem B(iii,b,i)"
print "Press enter to begin Problem B(iii,b,ii)"
my.pause()

#calculate X_squared and P for consecutive bin placement
sXsq1 = Xsq(10, r1Same)
sXsq2 = Xsq(10, r2Same)
sXsq3 = Xsq(10, r3Same)
sXsq4 = Xsq(10, j1Same)
sXsq5 = Xsq(10, j2Same)
sXsq6 = Xsq(10, j3Same)
sp1 = pval(10, sXsq1)
sp2 = pval(10, sXsq2)
sp3 = pval(10, sXsq3)
sp4 = pval(10, sXsq4)
sp5 = pval(10, sXsq5)
sp6 = pval(10, sXsq6)

sameXsq = [sXsq1, sXsq2, sXsq3, sXsq4, sXsq5, sXsq6]
sameP = [sp1, sp2, sp3, sp4, sp5, sp6]

#calculate X_squared and P for neighboring bin placement
aXsq1 = Xsq(10, r1Next)
aXsq2 = Xsq(10, r2Next)
aXsq3 = Xsq(10, r3Next)
aXsq4 = Xsq(10, j1Next)
aXsq5 = Xsq(10, j2Next)
aXsq6 = Xsq(10, j3Next)
ap1 = pval(10, aXsq1)
ap2 = pval(10, aXsq2)
ap3 = pval(10, aXsq3)
ap4 = pval(10, aXsq4)
ap5 = pval(10, aXsq5)
ap6 = pval(10, aXsq6)

adjaXsq = [aXsq1, aXsq2, aXsq3, aXsq4, aXsq5, aXsq6]
adjaP = [ap1, ap2, ap3, ap4, ap5, ap6]

#print results
print 'X_squared for consecutive bin placement:\n', sameXsq, '\n'
print 'Probability to exceed for consecutive bin placement:\n', sameP, '\n'
print 'X_squared for neighboring bin placement:\n', adjaXsq, '\n'
print 'Probability to exceed for neighboring bin placement:\n', adjaP, '\n'

print "End of Problem B(iii,b,ii)"
print "Press enter to end project"
my.pause()
