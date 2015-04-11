
# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab inline')


# In[2]:

from qutip import *
from scipy import *
from numpy import *
from scipy.optimize import minimize
import random
from math import log

from IPython.display import display

from sympy.interactive import printing
printing.init_printing()

from collections import defaultdict

from __future__ import division

from pygame import mixer
mixer.init() #you must initialize the mixer

from random import randrange
from datetime import datetime

from tempfile import TemporaryFile
arglikelihood = TemporaryFile()
argFid = TemporaryFile()


# In[3]:

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

pauliMatr = [[sz,sz],[sx,sx],[I,sx],[sz,I]]


######################################################
#GATE and DIMENSION
####################################################
N = 5        # n of qubits

G = toffoli()

interactions = [[0,1],[0,2],[1,2],[0,3],[1,3],[2,3],[0,4],[2,4],[1,4]]

#######################################################################
#USEFUL DEFINITIONS
#######################################################################

dimG = G.shape[0]  
CareStateDim = int(log(dimG,2))

h = N-CareStateDim
if h > 1 :
    dontCareStates = [sin(0.847)*basis(2,1) + cos(0.847)*basis(2,0)]*(h)
    dontCareIdentity = [qeye(2)]*(h)
else :
    dontCareStates = [sin(0.847)*basis(2,1) + cos(0.847)*basis(2,0)]
    dontCareIdentity = [qeye(2)]


dCS = tensor(dontCareStates)
dCI = tensor(dontCareIdentity)

###########################################################################
#FUNCTIONS
###########################################################################
from pygame import mixer
mixer.init() #you must initialize the mixer
alert = mixer.Sound('small-bell-ringing-02.wav') 

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
    
    i = 0
    k = 0
    H = 0
    
    for p in interactions :
        
        temp = 0
        
        for S in pauliMatr :
        
            OpChain = [qeye(2)]*N
            OpChain[p[0]] = S[0].op
            OpChain[p[1]] = S[1].op
            
            temp += x[k]*tensor(OpChain)
            k += 1

        H += temp 
        
    
    return H


# In[8]:

###########################################
##Functions to implement the Fidelity
###########################################

def getGate (G): #extract non zero elements of the gate and save them in s
    
    s = []

    rows = G.shape[0]
    colums = G.shape[1]
    
    for i in range(rows):
        for j in range(colums):
            if G[i][0][j] != 0 :
                s.append([i,j])
    return s

def get_bin(a):  #binary representation of a number a: useful to write the computational basis 
    
    s = bin(a)[2:]
    l = len(s)
    if l < 4 :
        r = '0'*(4-l)
        s = r + s
        
    return s

def getBasis (a) : #get the basis states according to the binary of a: 10 -> |10>
    
    if a == 0:
          B = [basis(2,0)]*(CareStateDim)
          return tensor(B)
    
    c = get_bin(a)
    if dimG == 16:
        return tensor(basis(2,int(c[0])) , basis(2,int(c[1])), basis(2,int(c[2])), basis(2,int(c[3])))
    if dimG == 8:
        return tensor(basis(2,int(c[0])) , basis(2,int(c[1])), basis(2,int(c[2])))
    if dimG == 4:
        return tensor(basis(2,int(c[1])) , basis(2,int(c[2])))
    
def Fidelity (J):  #Fidelity function
    
    s = getGate(G)
    H = HamiltonianAB(J)
    Fid = 1/(dimG + 1)
    U = (-j*H).expm()
    Udag = (j*H).expm()
    
    for x in s :
        for y in s:
            
            #definition of the basis kets and bras.             
            bra_i = getBasis(x[0]).dag()
            ket_j = getBasis(y[0])
            ket_k = getBasis(x[1])
            bra_l = getBasis(y[1]).dag()
            
            Epsilon = U*tensor(ket_k*bra_l, dCS*dCS.dag())*Udag
            Eps_ijkl = bra_i*(Epsilon.ptrace([0,1,2]))*ket_j
            
            Gstar_ik = G[x[0],x[1]].conjugate()
            G_jl = G[y[0],y[1]]
            
            
            fidStep = (1/(dimG*(dimG+1)))*Gstar_ik*Eps_ijkl*G_jl     
            Fid += abs(fidStep[0][0][0])
            
    return -Fid


# In[ ]:

########################
#SGD
######################

step = 0.0008

t = open('ToffoliOK', 'w+')
startTime = datetime.now()
Jopt = rand(len(interactions)*len(pauliMatr))

'''
Jopt = [-1.5891812232412459, 0.222201067110399, 2.9659134724454628, -0.61951776957560623, 
        0.61672394748136605, 2.9300251298898132, -2.8212390258083273, -0.87573248099454737, 
        0.40545591119003588, 0.30810380881645394, -0.20518203145596964, 2.1667996986171847, 
        0.012150621188495705, -0.78460485860462426, 0.48542689192214367, -0.2433259463869461, 
        -0.032646187693451088, 1.3827524967742615, -1.1725443013381454, 2.4916920168135519, 
        -0.00211546362361537, 0.56803067543107266, 0.15565154202542081, -0.95621213465733401, 
        -0.0081518495874134192, 0.70974885090792139, 0.27313441307128938, 1.3887347262985441, 
        -0.011186727441676907, 1.4691967774428702, 0.76662282735320775, 0.11736393526330254, 
        -0.041449262207567866, 0.10761097335491261, -0.051045946606989107, 3.2261405436903159]
'''


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
alert.play()


#---- NUN TE SCORDA' DE METTE ER MENO NELLA FIDELITY !!!!! ------#


from scipy.optimize import minimize
res = minimize(Fidelity, Jopt, method= "nelder-mead")


# In[7]:

Fidelity(rand(len(interactions)*len(pauliMatr)))


# In[ ]:




# In[9]:

i = s = 0
N = 5

for f in pauliMatr :
    s += 1
    c = 'J_' + f[0].name + f[1].name
    J = zeros(shape=(N,N))
    for k in interactions :
        J[k[0],k[1]] = Jopt[i]
        i += len(pauliMatr)
    i = s
    exec(c + '= asmatrix(J)')

