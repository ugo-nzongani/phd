import numpy as np
from qiskit import *
import random
from qiskit.quantum_info import Statevector

# https://arxiv.org/abs/2306.11747

def fidelity(psi1,psi2,n):
    F=0+0j
    for i in range(2**n):
        F+=np.conj(psi1[i])*psi2[i]
    return(abs(F)**2)
    
def compute_lk(m):
    l = []
    n = int(np.ceil(np.log2(m)))
    binary_m = bin(m)[2:].zfill(n)[::-1]
    k = 0
    for i in range(n):
        if binary_m[i] == '1':
            l.append(k)
        k += 1
    return l

def compute_M(i,l):
    if i == 0:
        return 2**l[i]
    else:
        return compute_M(i-1,l) + 2**l[i]
    
def list_M(l):
    res = []
    for i in range(len(l)-1):
        res.append(compute_M(i,l))
    return res

def uniform(m):
    if m == 1:
        print('m must be greater than 1')
    else:
        n = int(np.ceil(np.log2(m))) # number of qubits
        q = QuantumRegister(n,name='q')
        qc = QuantumCircuit(q)
        # true if m is a power of two
        if (m & (m-1) == 0) and m != 0:
            qc.h(q)
        else:
            l = compute_lk(m)
            k = len(l)-1
            m_list = list_M(l)

            for i in l[1:]:
                qc.x(q[i])

            m_0 = m_list[0]

            if l[0] > 0:
                for i in range(l[0]):
                    qc.h(q[i])
            theta_0 = -2*np.arccos(np.sqrt(m_0/m))
            qc.ry(theta_0,q[l[1]])

            i = l[0]
            while(i < l[1]):
                qc.ch(q[l[1]],q[i],ctrl_state='0')
                i += 1

            for i in range(1,k):
                theta_m = -2*np.arccos(np.sqrt(2**l[i]/(m-m_list[i-1])))
                qc.cry(theta_m,q[l[i]],q[l[i+1]],ctrl_state='0')
                j = l[i]
                while(j < l[i+1]):
                    qc.ch(q[l[i+1]],q[j],ctrl_state='0')
                    j += 1
        return qc
    
def main():
    m = random.randint(2, 10000000) # number of basis states
    qc = uniform(m)
    #qc.draw('mpl')
    
    state = Statevector.from_instruction(qc).data
    n = int(np.ceil(np.log2(m)))
    # target state
    target = np.zeros(2**n,dtype=np.complex128)
    for i in range(m):
        target[i] = 1/np.sqrt(m)

    # fidelity between the prepared and the target states
    print('Generation of '+str(m)+' basis states with '+str(n)+' qubits')
    print('Fidelity: ',fidelity(state,target,n))
    
if __name__ == "__main__":
    main()
