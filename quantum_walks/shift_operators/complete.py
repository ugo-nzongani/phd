from qiskit.circuit import *
from qiskit.circuit.library import MCXGate

# https://arxiv.org/abs/2304.01501

def shift_complete_graph(n):
    """Generates the quantum circuit implementing the shift operator of the complete graph K_{2**n}
    ----------
    n : int
        Number of qubits encoding the position

    Returns
    -------
    qiskit.circuit.quantumcircuit.QuantumCircuit
    """
    q = QuantumRegister(n,name='q')
    c = QuantumRegister(n,name='c')
    qc = QuantumCircuit(q,c)
    for i in range(n):
        qc.cx(c[i],q[i])
    return qc
