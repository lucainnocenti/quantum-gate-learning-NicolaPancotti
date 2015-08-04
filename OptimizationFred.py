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
N = 4        # n of qubits

GateMatrix = matrix([[1,0,0,0,0,0,0,0],
                     [0,1,0,0,0,0,0,0],
                     [0,0,1,0,0,0,0,0],
                     [0,0,0,1,0,0,0,0],
                     [0,0,0,0,1,0,0,0],
                     [0,0,0,0,0,1,0,0],
                     [0,0,0,0,0,0,1,0],
                     [0,0,0,0,0,0,0,-1]]) #GATE ODD
#G = Qobj(GateMatrix, dims = [[2,2,2],[2,2,2]])
#G = toffoli()
G = fredkin()
#######################################################################
#USEFUL DEFINITIONS
#######################################################################

dimG = G.shape[0]  
sq1 = float(sys.argv[2])
sq2 = float(sys.argv[3])

dCS = sqrt(sq1)*basis(2,1) + sqrt(sq2)*basis(2,0)
dCS = tensor(dCS)

########################
#SGD
######################

step = float(sys.argv[1])/10000

t = open('fredkin'+str(step)+'coef'+str(sq1)+str(sq2), 'w+')
startTime = datetime.now()

Jopt = rand(14)


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
        
        rho_0 = rand_ket(N = 8, dims = [[2,2,2], [1,1,1]])
        index = randrange(len(J))
        
        JdJ = [x for x in J]
        JdJ[index] += delta
        HJ = thesisLib.Hfred(J,N)
        HJdJ = thesisLib.Hfred(JdJ,N)
        grad = (thesisLib.Likelihood3(JdJ, rho_0,G,dCS,HJdJ) - thesisLib.Likelihood3(J, rho_0,G,dCS,HJ))/delta
        J[index] = J[index] + grad/sqrt(walk)
        HJ = thesisLib.Hfred(J,N)
        s = thesisLib.Likelihood3(J,rho_0,G,dCS,HJ) 
        t.write(str(time/step)+'    '+str(s)+'\n')        
        
        if abs(s) > 0.99 :
            check += 1    
            if check == 20:
                break
         
        if time/step % 100 < 1:
            print(str(time/step)+ '   ' + str(s)+ '   ' + str(walk))
            
            
        if time/step > 1600*(i+1) :
            print('HIT')
            break
    
    Jopt = J


print(datetime.now()-startTime)
print(Jopt)
t.close()
