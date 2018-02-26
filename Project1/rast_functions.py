import time
import numpy as np
#
#----------------------------------------------------------------------------
#
def loadrandom1a(num,high):
#
#define fixed size array and fill it with numbers input from keyboard
#num=total number of entries to be read from command line
#high=maximum value of positive integer entry
#
    random1=np.zeros(num,np.float64)
    print 'Enter '+str(num)+' integers on interval [0,'+str(high)+']:'
    for n in range(0,num):
        numin=input()
        if (numin < 0) or (numin > high):
            invalid_input=True
            print 'Integer must lie in interval [0,'+str(high)+']:'
            while invalid_input:
                numin=input()
                if (numin >= 0) and (numin <= high):
                    invalid_input=False
                else:
                    print 'Integer must lie in interval [0,'+str(high)+']:'
        random1[n]=numin
        print random1[0:num]
#
#    print random1.dtype
    return random1
#
#----------------------------------------------------------------------------
#
def loadrandom1b(num,high):
#
#define dynamic array and fill it
#num=total number of entries to be read from command line
#high=maximum value of positive integer entry
#
    random1=[]
    print 'Enter '+str(num)+' integers on interval [0,'+str(high)+']:'
    for n in range(0,num):
        numin=np.float64(input())
        if (numin < 0) or (numin > high):
            invalid_input=True
            print 'Integer must lie in interval [0,'+str(high)+']:'
            while invalid_input:
                numin=input()
                if (numin >= 0) and (numin <= high):
                    invalid_input=False
                else:
                    print 'Integer must lie in interval [0,'+str(high)+']:'
        random1.append(numin)
        print random1[0:num]
#
# Convert to a numpy arrray rather than a list (for performance, compactness, and functionality)
    random1=np.array(random1,np.float64)
#    print random1.dtype
    return random1
#
#----------------------------------------------------------------------------
#
def loadrandom1(num,high):
#
# Load numbers into an array from the keyboard.
    random1=np.zeros(num,np.float64)
    print 'Enter '+str(num)+' integers on interval [0,'+str(high)+']:'
    for n in range(0,num):
        numin=input()
        if (numin < 0) or (numin > high):
            invalid_input=True
            print 'Integer must lie in interval [0,'+str(high)+']:'
            while invalid_input:
                numin=input()
                if (numin >= 0) and (numin <= high):
                    invalid_input=False
                else:
                    print 'Integer must lie in interval [0,'+str(high)+']:'
        random1[n]=numin
#
    return random1
#
#----------------------------------------------------------------------------
#
def sleeper():
#
# introduces a user specified wait time
#
     print '\r'
     num = input('How long to wait: ')
     print '\r'
     time.sleep(num)
#
#----------------------------------------------------------------------------
#
def pause():
#
# pauses until a keyboard entry (e.g. carrage return)
#
     print '\r'
     dummy = raw_input('Pause')
     print '\r'
#
#----------------------------------------------------------------------------
#
def save(filename,varnames,vars):
    """
    filename: name of file, as a string
    varnames: name of variables to save, as list of strings
    vars: arrays of data to write to file, as list of data arrays
    USAGE: save('example.dat',['a','b','c'],[a,b,c])
    """
#
    file=open(filename,'wb')
#
    for i, name in enumerate(varnames):
        file.write(name+'\n')
        var=vars[i]
        shape=var.shape
        shape=','.join(np.array(shape,dtype=str))
        file.write(shape+'\n')
        type=str(var.dtype)
        file.write(type+'\n')
        svar=var.flatten().tobytes() #tobytes creates string representing var as stored in memory
        file.write(svar+'\n\n')
#
    file.close()
#
#----------------------------------------------------------------------------
#
def restore(filename):
    """
    filename: name of file to restore, as a string
    USAGE: x,y,... = restore('example.dat')
    """
#
    file=open(filename,'rb')
    data=[]
#
    print('Restoring variables')
    while True:
        varname=file.readline().replace('\n','')
        if varname == '': break        #escape if empty string is read (end of file)
        shape=file.readline().replace('\n','')
        shape=tuple(np.array(shape.split(','),dtype=int))
        type=file.readline().replace('\n','')
        svar=''
        line=''
        while line != '\n':
            svar += line
            line=file.readline()
        var=np.fromstring(svar[0:-1],dtype=type)
        var=var.reshape(shape)
        data.append(var)
#
    file.close()
    data=np.squeeze(data)       #remove leading dimension of 1 when restoreing a single array
    return data
#
#----------------------------------------------------------------------------
#
