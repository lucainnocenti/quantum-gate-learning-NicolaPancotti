from qutip import *
def Likelihood2(J, rho_0, G, dCS, H) : #Likelihood function
    
    j = complex(0,1)
    Ak = G*rho_0
    A = Ak*Ak.dag()
    rho = tensor(rho_0, dCS)
    U = (-j*H).expm()
    Btemp = U*rho
    Bk = Btemp*Btemp.dag()
    B = Bk.ptrace([0,1])
    out = (A*B).tr()

    return abs(out)

def Likelihood3(J, rho_0, G, dCS, H) : #Likelihood function
    j = complex(0,1)
    Ak = G*rho_0
    A = Ak*Ak.dag()
    rho = tensor(rho_0, dCS)
    U = (-j*H).expm()
    Btemp = U*rho
    Bk = Btemp*Btemp.dag()
    B = Bk.ptrace([0,1,2])
    out = (A*B).tr()

    return abs(out)



def Hfred(x,N) :
    k = 0
    H = 0
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()
    for q in [sx,sy,sz]:    
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[2] = q
        OpChain[1] = q
        temp += x[k]*tensor(OpChain)
        k+=1        
        H += temp 
    for p in [0,3]:
        for q in [sx,sz]:
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[2] = q
            OpChain[p] = q
            temp += x[k]*tensor(OpChain)
            OpChain = [qeye(2)]*N
            OpChain[1] = q
            OpChain[p] = q
            temp += x[k]*tensor(OpChain)
            k += 1
            H += temp 
    for q in [sx,sz]:
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[0] = q
        OpChain[3] = q
        temp += x[k]*tensor(OpChain)
        k+=1        
        H += temp 
    for i in range(4) :
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[i] = sz
        temp += x[k]*tensor(OpChain)
        H += temp 
        k += 1
    
    temp = 0
    OpChain = [qeye(2)]*N
    OpChain[3] = sx
    temp += x[k]*tensor(OpChain)#last one

    H += temp 

    
    return H


def H(x,N) :
    k = 0
    H = 0
    sx = sigmax()
    sz = sigmaz()
    for q in [sx,sz]:    
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[0] = q
        OpChain[1] = q
        temp += x[k]*tensor(OpChain)
        k+=1        
        H += temp 
    for p in [2,3]:
        for q in [sx,sz]:
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[0] = q
            OpChain[p] = q
            temp += x[k]*tensor(OpChain)
            OpChain = [qeye(2)]*N
            OpChain[1] = q
            OpChain[p] = q
            temp += x[k]*tensor(OpChain)
            k += 1
            H += temp 
    for q in [sx,sz]:
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[2] = q
        OpChain[3] = q
        temp += x[k]*tensor(OpChain)
        k+=1        
        H += temp 
    for i in range(4) :
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[i] = sz
        temp += x[k]*tensor(OpChain)
        H += temp 
        k += 1
    
    temp = 0
    OpChain = [qeye(2)]*N
    OpChain[3] = sx
    temp += x[k]*tensor(OpChain)#last one
    H += temp 

    
    return H

def HchainHeis(x,N) :
    
    k = 0
    H = 0

    sx = sigmax()
    sz = sigmaz()
    sy = sigmay()

    for p in range(1,N-1) :
        for q in [sx,sy,sz] :
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[p] = q
            OpChain[p+1] = q
            
            temp += x[k]*tensor(OpChain)
            k += 1
            H += temp 
    for q in [sx,sy,sz] :
        temp = 0
        OpChain = [qeye(2)]*N
        OpChain[0] = q
        OpChain[N-1] = q
        
        temp += x[k]*tensor(OpChain)
        k += 1
        H += temp 
    for i in range(N):
        for q in [sx,sy,sz] :
            temp = 0
            OpChain = [qeye(2)]*N
            OpChain[i] = q
            temp += x[k]*tensor(OpChain)
            k += 1
            H += temp

            
    
    return H

