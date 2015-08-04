import sys
sys.path.append("/afs/ipp-garching.mpg.de/home/n/nicolap/lib/python/")
from qutip import *
from scipy import *
from numpy import *
import random
import thesisLib
from math import log
from random import randrange
from datetime import datetime

j = complex(0,1)

######################################################
#GATE and DIMENSION
####################################################
N = int(sys.argv[4])        # n of qubits
G = cnot()

#######################################################################
#USEFUL DEFINITIONS
#######################################################################

dimG = G.shape[0]  
CareStateDim = int(log(dimG,2))
h = N-CareStateDim
sq1 = float(sys.argv[2])
sq2 = float(sys.argv[3])

if h > 1 :
    dontCareStates = [sq1*basis(2,1) + sq2*basis(2,0)]*(h)
    dontCareIdentity = [qeye(2)]*(h)
else :
    dontCareStates = [sq1*basis(2,1) + sq2*basis(2,0)]
    dontCareIdentity = [qeye(2)]


dCS = tensor(dontCareStates)
dCI = tensor(dontCareIdentity)

########################
#SGD
######################

step = float(sys.argv[1])/10000

t = open('chainOut/cnotHeis'+str(step)+'coef'+str(sq1)+str(sq2)+'N'+str(N), 'w+')
startTime = datetime.now()

Jopt = rand(3*(2*N-1))


delta = 0.0001
fcheck = 0
s = 0
time = step


for i in range(3):
    
    J = [x for x in Jopt]
    
    walk = step*500*i
    
    while True :
        
        walk += step
        time += step
        
        rho_0 = rand_ket(N = 4, dims = [[2,2], [1,1]])
        index = randrange(len(J))
        
        JdJ = [x for x in J]
        JdJ[index] += delta
        HJ = thesisLib.HchainHeis(J,N)
        HJdJ = thesisLib.HchainHeis(JdJ,N)
        grad = (thesisLib.Likelihood2(JdJ, rho_0,G,dCS,HJdJ) - thesisLib.Likelihood2(J, rho_0,G,dCS,HJ))/delta
        J[index] = J[index] + grad/sqrt(walk)
        HJ = thesisLib.HchainHeis(J,N)
        s = thesisLib.Likelihood2(J,rho_0,G,dCS,HJ) 
        t.write(str(time/step)+'    '+str(s)+'\n')        
        
        if abs(s) > 0.99 :
            check += 1    
            if check == 20:
                break
         
        #if time/step % 100 < 1:
         #   print(str(time/step)+ '   ' + str(s)+ '   ' + str(walk))
            
            
        if time/step > 1600*(i+1) :
            #print('HIT')
            break
    
    Jopt = J


print(datetime.now()-startTime)
print(Jopt)
t.close()

