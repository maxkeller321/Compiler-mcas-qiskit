from qiskit import QuantumCircuit
from qiskit.compiler import transpile
from qiskit.compiler.transpile import CouplingMap
from header_footer_mcas import mcas_footer, mcas_header
import numpy as np 
import textwrap
import os 


mcas_file_name = "mcas_test.py"
indent_level = "\t\t\t"
            

def transpile_ciruit_for_diamond(quantum_circuit_instance):
    """
    Compiles/ transpiles the QuantumCircuit instance from qiskit to a quantum_circuit_instance instance which can be run on our NV-Center QC
        params: 
            quantum_circuit_instance: instance of quantum_circuit_instance
        return: 
            transpiled_qc: transpiled QunatumCircuit instance
    """

    if not isinstance(quantum_circuit_instance, QuantumCircuit): 
        raise Exception("The passed instance (quantum_circuit_instance) must come from the qiskit class QuantumCircuit")
    
    coupling_string = [[0, 1], [1, 0], [2, 0], [0,2], [1,2], [2, 1]]
    CM = CouplingMap(coupling_string)
    basis_gates = ['id', 'ry', 'rx', 'cx']
    transpiled_qc = transpile(quantum_circuit_instance, coupling_map=CM, basis_gates=basis_gates, optimization_level=3, seed_transpiler=1)

    return transpiled_qc


def construct_list_with_json_from_QuantumCircuit_instance(quantum_circuit_instance): 
    """
        constructs a list of json objects from a qiskit QuantumCircuit instance
        where each json object contains the important data for a single gate/ operation
        The order in the list is equivalent to the execution order for the quantum computer.
        params: 
            quantum_circuit_instance: instance of quantum_circuit_instance
        return: 
            list of json objects where each object contains the important data for a single gate/ operation
    """

    if not isinstance(quantum_circuit_instance, QuantumCircuit): 
        raise Exception("The passed instance (quantum_circuit_instance) must come from the qiskit class QuantumCircuit")
    
    quantum_circuit_operations = []
    quantum_circuit_data = quantum_circuit_instance.data

    for gate_data in quantum_circuit_data:
        operation_json = {}

        gate_instance = gate_data[0]
        operation_name = gate_instance.name 
        gate_params = gate_instance.params

        operation_json['operation_name'] = operation_name
        operation_json['operation_params'] = gate_params
        

        if gate_instance.num_qubits == 2: 
            controlled_qubit = gate_data[1][1].index
            controlling_qubit = gate_data[1][0].index
            operation_json['controlled_qubit'] = controlled_qubit
            operation_json['controlling_qubit'] = controlling_qubit 
        else: 
            qubit_index = gate_data[1][0].index
            operation_json['qubit_index'] = qubit_index

        
        quantum_circuit_operations.append(operation_json)

    return quantum_circuit_operations


class McasWriter:
    def __init__(self, gate_operation_json):
        self.mcas_line_list = [] # list of strings with the operations for the mcas file 
    
    def add_rx_nitrogen(self, theta, ms=0):
        self.mcas_line_list.append("rx_nitrogen(mcas, {}, ms={})".format(theta, ms))
    def add_ry_nitrogen(self, theta, ms=0):
        self.mcas_line_list.append("ry_nitrogen(mcas, {}, ms={})".format(theta, ms))
    def add_rx_c_414(self, theta, ms=0):
        self.mcas_line_list.append("rx_carbon_414(mcas, {}, ms={})".format(theta, ms))
    def add_ry_c_414(self, theta, ms=0):
        self.mcas_line_list.append("ry_carbon_414(mcas, {}, ms={})".format(theta, ms))
    def add_rx_c_90(self, theta, ms=0):
        self.mcas_line_list.append("rx_carbon_90(mcas, {}, ms={})".format(theta, ms))
    def add_ry_c_90(self, theta, ms=0):
        self.mcas_line_list.append("ry_carbon_90(mcas, {}, ms={})".format(theta, ms))
    def add_cx(self, controlled_qubit, controlling_qubit): 
        print("lalalal")
    def add_nuclear_spin_initialisation(self, state): 
        print("lalalla")
    def add_nuclear_spin_readout(self, state):
        print("lalalla")

    def write_mcas_file(self, destination): 
        with open(os.path.join(destination, mcas_file_name), 'w') as mcas_file: 
            
            mcas_file.write(mcas_header)

            # Operation content 
            for line in self.mcas_line_list:
                indented_line = textwrap.indent(text=line, prefix=indent_level) 
                mcas_file.writelines(indented_line + '\n')
            
            mcas_file.write(mcas_footer)

