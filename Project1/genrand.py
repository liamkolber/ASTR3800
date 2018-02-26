import random as rand


random_nums = []
ht = []

for i in range (10):
	random_nums.append(rand.randrange(1,10,1))
print (random_nums), '\n'


for i in range (10):
	ht.append(rand.randrange(1,3,1))
	if ht[i] == 1:
		ht[i] = "H"
	else:
		ht[i] = "T"
print (ht)
