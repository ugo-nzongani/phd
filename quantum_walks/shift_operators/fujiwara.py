from qiskit import *

# We use the convention that the walkers goes -1 if the coin is |1> and +1 if it is |0>

def shift_gates(n):
    """Quantum gates of Z+ and Z- used to move the walker
    Parameters
    ----------
    n : int
        The number of qubits encoding the position

    Returns
    -------
    qiskit.circuit.gate.Gate,qiskit.circuit.gate.Gate
    """
    # Position register
    b = QuantumRegister(n, name="b")
    # Coin register
    s = QuantumRegister(1, name="s")
    qc = QuantumCircuit(b, s)
    for i in range(n):
        ctrl = [s[0]] + b[:n-i-1]
        qc.mcx(ctrl, [b[n-i-1]])
    z_up = qc.to_gate(label="Z+")
    z_down = qc.inverse().to_gate(label="Z-")
    return z_up, z_down

def fujiwara_shift(n):
    """Quantum circuit implementing the shift operator introduced by Fujiwara et al. (10.1103/PhysRevA.72.032329)
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
    qc = QuantumCircuit(b, s)
    qubits = [i for i in b] + [s[0]]
    z_up, z_down = shift_gates(n)
    qc.x(s[0])
    qc.append(z_up, qubits)
    qc.x(s[0])
    qc.append(z_down, qubits)
    return qc
