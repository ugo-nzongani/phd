from qiskit.circuit import *
import numpy as np
from qiskit.circuit.library import MCXGate

# https://arxiv.org/abs/2304.01501

def generate_gray_code(n, binary=False):
    """Generates a gray list of 2**n integers
    Parameters
    ----------
    n : int
        Binary integer to convert
    binary : bool
        Set to False if the list coefficient should be integers instead of strings
        
    Returns
    -------
    str array or int array
    """
    if n < 0:
        raise ValueError("Input must be a non-negative integer.")

    if n == 0:
        return ['0'*n] if binary else [0]
    
    gray_list = ['0'*n, int_to_binary(1,n,reverse=False)] if binary else [0, 1]
    
    for i in range(1, n):
        mirror = 2**i
        if binary:
            gray_list += [bin(mirror + int(num, 2))[2:].zfill(n) for num in reversed(gray_list)]
        else:
            gray_list += [mirror + int(num) for num in reversed(gray_list)]
    
    return gray_list

def int_to_binary(k,n,reverse=True):
    """Converts an integer k to its binary representation on n bits
    Parameters
    ----------
    k : int
        Integer to convert
    n : int
        Number of bits on which k should be written
    reverse : bool
        If True the indices of the bits in the bit string of k are reversed, i.e. binary[0] corresponds to
        the least significant bit of k
        
    Returns
    -------
    str
    """
    binary = bin(k)[2:].zfill(n)
    if reverse:
        binary = binary[::-1]
    return binary
    
def get_coin_dim(n):
    """Returns the dimension of the coin
    ----------
    n : int
        Dimension of the hypercube

    Returns
    -------
    int
    """
    return int(np.ceil(np.log2(n)))

def shift_hypercube_power_of_two(n,coin_dim,gray_list):
    """Generates the quantum circuit implementing the shift operator of the n-hypercube with n a power
    of two
    ----------
    n : int
        Dimension of the hypercube
    coin_dim : int
        Dimension of the coin
    gray_list : int/str list
        Gray code list

    Returns
    -------
    qiskit.circuit.quantumcircuit.QuantumCircuit
    """
    q = QuantumRegister(n,name='q')
    c = QuantumRegister(coin_dim,name='c')
    qc = QuantumCircuit(q,c)
    for i,j in enumerate(gray_list):
        gate = MCXGate(coin_dim,ctrl_state=j)
        qc.append(gate,qc.qubits[n:]+qc.qubits[i:i+1])
    return qc

def shift_hypercube(n,coin_dim):
    """Generates the quantum circuit implementing the shift operator of the n-hypercube with self-loops
    if n is not a power of two
    ----------
    n : int
        Dimension of the hypercube

    Returns
    -------
    qiskit.circuit.quantumcircuit.QuantumCircuit
    """
    gray_list = generate_gray_code(coin_dim,binary=True)

    q = QuantumRegister(n,name='q')
    c = QuantumRegister(coin_dim,name='c')
    qc = QuantumCircuit(q,c)

    if (n & (n-1) == 0) and n != 0:
        qc.append(shift_hypercube_power_of_two(n,coin_dim,gray_list),qc.qubits)
    else:
        next_power_of_two = int(2**(coin_dim))
        k = next_power_of_two-n
        qc.append(shift_hypercube_power_of_two(n,coin_dim,gray_list[k:]),qc.qubits)

    return qc
