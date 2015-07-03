import sys
sys.path.append("/afs/ipp-garching.mpg.de/home/n/nicolap/lib/python/")
from qutip import *
from scipy import *
from numpy import *
import random
from math import log
from random import randrange
from datetime import datetime

j = complex(0,1)

class PAULI (object) :
    
    name = ''
    op = Qobj()
    
    def __init__(self, name, PauliMatrix) :
        self.name = name
        self.op = PauliMatrix
    
def make_pauli_matrix (name, PauliMatrix) :
    Pmatrix = PAULI (name, PauliMatrix)
    return Pmatrix

 
sx = make_pauli_matrix('X', sigmax())
sy = make_pauli_matrix('Y', sigmay())
sz = make_pauli_matrix('Z', sigmaz())
I = make_pauli_matrix('I', qeye(2))

######################################################
#GATE and DIMENSION
####################################################
N = 5        # n of qubits

GateMatrix = matrix([[1,0,0,0,0,0,0,0],
                     [0,1,0,0,0,0,0,0],
                     [0,0,1,0,0,0,0,0],
                     [0,0,0,1,0,0,0,0],
                     [0,0,0,0,1,0,0,0],
                     [0,0,0,0,0,1,0,0],
                     [0,0,0,0,0,0,1,0],
                     [0,0,0,0,0,0,0,-1]]) #GATE ODD
#G = Qobj(GateMatrix, dims = [[2,2,2],[2,2,2]])
G = toffoli()

interactions = [[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4]]
particles = [0,1,2,3,4]

#######################################################################
#USEFUL DEFINITIONS
#######################################################################

dimG = G.shape[0]  
CareStateDim = int(log(dimG,2))

h = N-CareStateDim

sq1 = float(sys.argv[2])
sq2 = float(sys.argv[3])


if h > 1 :
    dontCareStates = [sqrt(sq1)*basis(2,1) + sqrt(sq2)*basis(2,0)]*(h)
    dontCareIdentity = [qeye(2)]*(h)
else :
    dontCareStates = [sqrt(sq1)*basis(2,1) + sqrt(sq2)*basis(2,0)]
    dontCareIdentity = [qeye(2)]


dCS = tensor(dontCareStates)
dCI = tensor(dontCareIdentity)

###########################################################################
#FUNCTIONS
###########################################################################
def Likelihood(J, rho_0) : #Likelihood function
    
    Ak = G*rho_0
    #A = tensor(Ak*Ak.dag(), dCI)
    A = Ak*Ak.dag()
    
    rho = tensor(rho_0, dCS)
    H = HamiltonianAB(J)
    U = (-j*H).expm()
    Btemp = U*rho
    Bk = Btemp*Btemp.dag()
    B = Bk.ptrace([0,1,2])
    
    out = (A*B).tr()

    return abs(out)

def HamiltonianAB(x) :
    
    k = 0
    H = 0


    for q in [sx,sz]:    
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[0] = q.op
        OpChain[1] = q.op
        temp += x[k]*tensor(OpChain)
        k+=1        
        H += temp 

    for p in [2,3,4]:
        for q in [sx,sz]:
            
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[0] = q.op
            OpChain[p] = q.op
            
            temp += x[k]*tensor(OpChain)
            
            OpChain = [qeye(2)]*N
            OpChain[1] = q.op
            OpChain[p] = q.op
            
            temp += x[k]*tensor(OpChain)
            k += 1
            
            H += temp 
                
    for p in [3,4]:
        for q in [sx,sz]:

    
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[2] = q.op
            OpChain[p] = q.op
            temp += x[k]*tensor(OpChain)
            k+=1        
            H += temp 
    

    for i in range(2) :
    
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[i] = sz.op
        temp += x[k]*tensor(OpChain)
    
        H += temp 

    k += 1
    temp = 0
    OpChain = [qeye(2)]*N
    OpChain[2] = sx.op
    temp += x[k]*tensor(OpChain)#last one

    H += temp 

    
    return H



def Hamiltonian(x) :
    
    k = 0
    H = 0
    i = 0
    
    for p in interactions :
        for q in [sx,sz] :
        
            temp = 0
            
            OpChain = [qeye(2)]*N
            OpChain[p[0]] = q.op
            OpChain[p[1]] = q.op
            
            temp += x[k]*tensor(OpChain)
            k += 1
    
            H += temp 
    
    for i in range(2) :
    
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[i] = sz.op
        temp += x[k]*tensor(OpChain)
        k += 1
    
        H += temp 

    temp = 0
    OpChain = [qeye(2)]*N
    OpChain[2] = sx.op
    temp += x[k]*tensor(OpChain)

    H += temp 

    
    return H

########################
#SGD
######################

step = float(sys.argv[1])/10000

t = open('toffoli_optimized'+str(step)+'coef'+str(sq1)+str(sq2), 'w+')
startTime = datetime.now()

Jopt = rand(14)


delta = 0.0001
check = 0
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
        grad = (Likelihood(JdJ, rho_0) - Likelihood(J, rho_0))/delta
        J[index] = J[index] + grad/sqrt(walk)
        s = Likelihood(J, rho_0) 
        
        
        
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

