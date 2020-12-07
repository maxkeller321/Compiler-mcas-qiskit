from qiskit import QuantumCircuit
from qiskit.compiler import transpile
from qiskit.compiler.transpile import CouplingMap
import numpy as np 
import textwrap
import os 
from collections.abc import Iterable
from .header_footer_mcas import mcas_footer, mcas_header


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

def match_integers_with_nuclear_spin(qubit_integer):
    if qubit_integer == 0:
        nuclear_spin = '14n'
    elif qubit_integer == 1: 
        nuclear_spin = '13c414'
    elif qubit_integer == 2: 
        nuclear_spin = '13c90'
    else: 
        raise Exception("qubit_integer must be between 0<=x<=2, but was {}".format(qubit_integer))
    
    return nuclear_spin

class McasWriter:
    def __init__(self):
        self.mcas_line_list = [] # list of strings with the operations for the mcas file 

    def interprete_json_operations(self, gate_operation_json):
        if isinstance(gate_operation_json, Iterable):
            
            for operation in gate_operation_json:
                if 'qubit_index' in operation and 'operation_name' in operation and 'operation_params' in operation: 
                    qubit_index = operation['qubit_index'] 

                    if operation['operation_name'] not in ['barrier', 'measure']:
                        if qubit_index == 0: 
                            if operation['operation_name'] == 'rx': 
                                self.add_rx_nitrogen(operation['operation_params'][0])
                            elif operation['operation_name'] == 'ry': 
                                self.add_ry_nitrogen(operation['operation_params'][0])
                            else: 
                                raise Exception("Invalid operation encountered in qubit {} the operation json" \
                                            "Only rx, ry and cnot are valid for our Backend.".format(qubit_index))
                        elif qubit_index == 1: 
                            if operation['operation_name'] == 'rx': 
                                self.add_rx_c_414(operation['operation_params'][0])
                            elif operation['operation_name'] == 'ry': 
                                self.add_ry_c_414(operation['operation_params'][0])
                            else: 
                                raise Exception("Invalid operation encountered in qubit {} the operation json" \
                                            "Only rx, ry and cnot are valid for our Backend.".format(qubit_index))

                        elif qubit_index == 2: 
                            if operation['operation_name'] == 'rx': 
                                self.add_rx_c_90(operation['operation_params'][0])
                            elif operation['operation_name'] == 'ry': 
                                self.add_ry_c_90(operation['operation_params'][0])
                            else: 
                                raise Exception("Invalid operation encountered in qubit {} the operation json" \
                                            "Only rx, ry and cnot are valid for our Backend.".format(qubit_index))
                        else: 
                            raise Exception("Invalid qubit index encountered. Our system does only have 3 qubits: 0, 1 & 2")
            
                elif 'controlled_qubit' in operation and 'controlling_qubit' in operation and 'operation_name' in operation:
                    self.add_cx(operation['controlled_qubit'], operation['controlling_qubit'])

                elif 'operation_name' in operation: 
                    if operation['operation_name'] == 'measure': 
                        pass 
                    elif operation['operation_name'] == 'barrier': 
                        pass
                    else: 
                        raise Exception("The passed operation is not valid!")
                    
                else: 
                    raise Exception("The passed operation is not valid!")
        else: 
            raise Exception("{} is not iterable!. Insert a valid list with json operations!".format(gate_operation_json))
    
    def add_rx_nitrogen(self, theta, ms=0):
        self.mcas_line_list.append("rx_nitrogen(mcas, {}, ms={})".format(theta, ms))
    def add_ry_nitrogen(self, theta, ms=0):
        self.mcas_line_list.append("ry_nitrogen(mcas, {}, ms={})".format(theta, ms))
    def add_rx_c_414(self, theta, ms=-1):
        self.mcas_line_list.append("rx_carbon_414(mcas, {}, ms={})".format(theta, ms))
    def add_ry_c_414(self, theta, ms=-1):
        self.mcas_line_list.append("ry_carbon_414(mcas, {}, ms={})".format(theta, ms))
    def add_rx_c_90(self, theta, ms=-1):
        self.mcas_line_list.append("rx_carbon_90(mcas, {}, ms={})".format(theta, ms))
    def add_ry_c_90(self, theta, ms=-1):
        self.mcas_line_list.append("ry_carbon_90(mcas, {}, ms={})".format(theta, ms))
    def add_cx(self, controlled_qubit, controlling_qubit): 
        controlled_nuclear_spin = match_integers_with_nuclear_spin(controlled_qubit)
        controlling_nuclear_spin = match_integers_with_nuclear_spin(controlling_qubit)
        self.mcas_line_list.append("cnot_between_nuclear_spins(mcas, '{}', '{}')".format(controlled_nuclear_spin, controlling_nuclear_spin))
    def add_nuclear_spin_initialisation(self, state): 
        self.mcas_line_list.append("initialise_nuclear_spin_register(mcas, '{}')".format(state))
    def add_nuclear_spin_readout(self, state):
        self.mcas_line_list.append("readout_nuclear_spin_state(mcas, '{}')".format(state))

    def write_mcas_file(self, destination): 
        with open(os.path.join(destination, mcas_file_name), 'w') as mcas_file: 
            
            mcas_file.write(mcas_header)

            # Operation content 
            for line in self.mcas_line_list:
                indented_line = textwrap.indent(text=line, prefix=indent_level) 
                mcas_file.writelines(indented_line + '\n')
            
            mcas_file.write(mcas_footer)

