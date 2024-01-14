from qiskit import *
from qiskit.circuit.library import QFT
from numpy import pi

# We use the convention that the walker goes -1 if the coin is |1> and +1 if it is |0>

def shakeel_shift(n):
    """Shift operator for cycle of size power of two using the QFT: https://arxiv.org/abs/1912.00978
    Parameters
    ----------
    n : int
        The number of qubits encoding the position

    Returns
    -------
    qiskit.circuit.quantumcircuit.QuantumCircuit
    """
    # Position register
    b = QuantumRegister(n, name="b")
    # Coin register
    s = QuantumRegister(1, name="s")
    # Create a quantum circuit using the quantum and classical registers
    qc = QuantumCircuit(b, s)
    position_qubits = qc.qubits[:n]
    qft = QFT(n, do_swaps=False)
    #qc.x(s[0]) # uncomment if the walker goes +1 when the coin is |1>
    for i in range(n):
        qc.cx(s[0], b[n-i-1])
    qc.append(qft, position_qubits)
    for i in range(n):
        qc.p(2*pi/2**(n-i), n-i-1)
    qc.append(qft.inverse(), position_qubits)
    for i in range(n):
        qc.cx(s[0], b[i])
    #qc.x(s[0])
    return qc
