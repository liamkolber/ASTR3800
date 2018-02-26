import numpy as np
import rast_functions as my
#
print 'cointoss.py'
print 'Generating random numebers from specified seed:'
my.pause()
#
numelements=10
#
#Hard code seed.
#
seed=62659 
coin=np.chararray(numelements) #set array size
coin[:]='T' #initialize the array
#
#The next five lines of code yield one numelements long string of heads or tails
#
np.random.seed(seed)
rand=np.random.uniform(size=(numelements))
coin[np.where(rand < 0.5)]='H'
print 'seed: ',seed
print 'result: ',coin
#
my.pause()
#
print 'Input seed from keyboard:'
my.pause()
#
#Hard codeing the seed is good for code reproducability, but is inflexible. Instead, perhaps read the seed for the random number generator from the terminal.
#
print 'Seed?'
seed=input()
np.random.seed(seed)
rand=np.random.uniform(size=(numelements))
coin[np.where(rand < 0.5)]='H'
print 'seed: ',seed
print 'result: ',coin
#
my.pause()
#
print 'Iterate for a specific random sequnce:'
my.pause()
#
#The code below iterates many random number sequences until arriving at the specified H,T,H,T,etc. sequence, and counts how many trials it took.
#
n=0
seed=62659 
np.random.seed(seed)
#
#The two while statements below perform identically and the syntax choice is stylistic.
#while (coin != ['H','T','H','T','H','T','H','T','H','T']).any():
while any(coin != ['H','T','H','T','H','T','H','T','H','T']):
    coin[:]='T'
    state=np.random.get_state()
    rand=np.random.uniform(size=(numelements))
    coin[np.where(rand < 0.5)]='H'
    n=n+1
#    print n
#    print state
    print 'result: ',n,' ',coin
#
my.pause()
#
print 'Generate the same random sequence starting from a saved state of the random number generator:'
# The code below generated the same sequence by starting from the state of random  number generator
#
np.random.set_state(state)
rand=np.random.uniform(size=(numelements))
coin[np.where(rand < 0.5)]='H'
print 'result: ',coin
#
print '\r'
print 'cointoss.py'
