import numpy as np

#### Hamiltonian parsing and formatting ####
def parse_hamiltonian(filename, filepath):
    '''
    Takes a txt file with the Hamiltonian, where each line contains a linear coefficient (h_pq or h_pqrs) and Pauli string.
    Returns the coefficient ("scalar") and the Pauli string ("gates"), separately (in a tuple)
    In principle, this can be applied to any qubit operator, not just the Hamiltonian.
    '''
    parsed_data = []
    with open(filepath+'{}'.format(filename)) as file:
        for line in file:
            split_line = line.split(' ', 1)
            scalar = split_line[0]
            gates = split_line[1].split(']')[0] + ']'
            # print(scalar, gates)
            parsed_data.append((scalar, gates))
    return parsed_data

def format_hamiltonian(parsed_data):
    '''
    Formats the Hamiltonian or other qubit operator into a string that can be read by QuestLink.
    Careful: This is very manual, hard-coded, imperfect code. All that the code does is to read a string and transform it into another string.
    It might easily be that if the structure of the operator that you're parsing is
    different from the one that I used as an example to write this code, things might go wrong.
    Therefore, always check the formatted operator for bugs, or feel free to make the code more abstract so that it generalizes better. 
     '''
    formatted_data = '{' # one string, where each line is one Hamiltonian term
    
    for i in range(len(parsed_data)):
        scalar = parsed_data[i][0]
        gates = parsed_data[i][1]
        line = '{ '
        
        #### format the coefficient of the Hamiltonian ("scalar") ####
        if 'e' in scalar: # Mathematica doens't like the "e-10" notation for decimals
            e_letter_index = scalar.index('e')
            e_minus = scalar[e_letter_index:e_letter_index+2] # should keep the sign "e-" and replace with "10^-"
            scalar = scalar.replace(e_minus, '*10^-')
            
        line += scalar
        
        #### format the qubit operators in the Hamiltonian ("gates") ####
        if gates != '[]': # only proceed if there is actually a gate; for the constant hamiltonian term with no gates, do nothing
            newgates = gates.replace('[', '')
            newgates = newgates.replace(']', '')
            newgates = newgates.split(' ')

            for i in range(len(newgates)):
                pauli_letter = newgates[i][0] # the letter X, Y, or Z
                qubit_index = newgates[i][1:] # the index of the qubit on which the gate operates
                line += ', Subscript[' + pauli_letter + ', ' + qubit_index + ']'            
        line += '},'
        formatted_data += '\n' + line
    formatted_data = formatted_data[:-1] # this just removes the last comma
    formatted_data += '\n}'
    return formatted_data

def get_formatted_hamiltonian(filename, filepath):
    parsed_data = parse_hamiltonian(filename, filepath)
    formatted_ham = format_hamiltonian(parsed_data)
    return formatted_ham

#### Ansatz circuit parsing and formatting ####
def parse_circuit(filename, filepath):
    parsed_data = []
    with open(filepath+'{}'.format(filename)) as file:
        for line in file:
            split_line = line.rsplit(' ', 1)
            qubit_index_raw = split_line[1]
            qubit_index = ''.join(i for i in qubit_index_raw if i.isdigit()) # this is either one number (for single-qubit gates) or two numbers (2-qubit gates)
            gate = split_line[0]
            parsed_data.append((qubit_index, gate))
        print('Number of qubits:', parsed_data[0][0])
        del parsed_data[0] # the first line just defines a quantum register with a certain number of qubits (we just print that number here)
    return parsed_data

def format_circuit(parsed_data):
    '''
    This function takes the parsed raw circuit desinged with qiskit syntax and formats it into the QuestLink syntax.
    Again, it is very hard-coded and tailored to work for the kUpCCGSD ansatz, which has U1, U3 and CNOT gates as defined by qiskit.
    The code will not work for any arbitrary ansatz, although it could provide a starting point for a more general code.
    '''
    mapping = {'t0': 't0',
            't1': 't1',
            't2': 't2',
            't3': 't3',
            't4': 't4',
            't5': 't5',
            't6': 't6',
            't7': 't7',
            't8': 't8',
            't9': 't9',
            't10': 't10',
            't11': 't11',
            't0_1': 't12',
            't1_1': 't13',
            't2_1': 't14',
            't3_1': 't15',
            't4_1': 't16',
            't5_1': 't17',
            't6_1': 't18',
            't7_1': 't19',
            't8_1': 't20',
            't9_1': 't21',
            't10_1': 't22',
            't11_1': 't23',
            }

    # Should generalize this mapping function. It was written for the kUp ansatz

    formatted_circuit = '{\n' # a list of strings, where each element is one Hamiltonian term 

    for i in range(len(parsed_data)):
        qubit_index = parsed_data[i][0]
        gate = parsed_data[i][1]

        formatted_circuit += '{' # each element comes inside a curly bracket     
        if 'u' in gate:
            assert len(qubit_index) == 1 # single-qubit gate
            param = gate[3:] # take only the part after open paranthesis
            
            k = param.rfind(")")
            param = param[:k] + "" + param[k+1:] #  this just removes the final close paranthesis ")" 

            param = param.replace('pi', 'Pi')
            print('old: ', param)
            
            if 't' in param: # selects only the parametrized gates
                t_index = param.index('t')
                symbol_and_morecrap = param[t_index:t_index+5]
                symbol_only = ''
                for i in symbol_and_morecrap:
                    if i.isdigit() or i == 't' or i == '_':
                        symbol_only += i # symbol only really is the symbol t0, t1, ..., t1_23 and nothing else
                symbol_new = mapping[symbol_only] # t0, t1, ..., t12, t13, ..., t23
                number = symbol_new[1:] # get the number of the symbol only 0, 1, ..., 12, 13, 23
                symbol_questformat = 'Subscript[\[Theta]' + ', ' + number + ']' # this line is wrong if param is (theta - 1/2) * pi !
                param = param.replace(symbol_only, symbol_questformat)
                print('new: ', param)

            if 'u1' in gate:
                gatelabel = 'U1'
            elif 'u3' in gate:
                gatelabel = 'U3'
            term = 'Subscript[' + gatelabel + ', ' + qubit_index + ']' 
            term += ', ' + param
        elif 'cx' in gate:
            assert len(qubit_index) == 2 # two-qubit gate
            term = 'Subscript[' + 'C' + ', ' + qubit_index[0] + ']' + '[' + 'Subscript[' + 'X' + ', ' + qubit_index[1] + ']]'
        formatted_circuit += term
        formatted_circuit += '},\n'
    formatted_circuit = formatted_circuit[:-2] # this just removes the last comma and the \n
    formatted_circuit += '\n}'
    return formatted_circuit


def get_formatted_circuit(filename, filepath):
    parsed_data = parse_circuit(filename, filepath)
    formatted_circuit = format_circuit(parsed_data)
    return formatted_circuit

#### Write to textfile ####
def write_to_textfile(formatted_object, filename, path):
    filename = filename + '_formatted.txt'
    text_file = open(path+'{}'.format(filename), "w")
    n = text_file.write(formatted_object)
    text_file.close()


#### Examples ####

## Hamiltonian
#n_points = 40 
#bond_length_interval = 4.0 / n_points
#
#for point in range(1, 41):
#    bond_length = np.round(bond_length_interval * float(point) + 0.2, 2)
#    print(bond_length)
#    filename = 'ham_LiH_rhforbs_' + str(bond_length)
#    filepath = 'hamiltonians/LiH/' # this folder contains the 40 txt files with the LiH Hamiltonian in its raw (qiskit) format
#    writepathname = 'hamiltonians/LiH/formatted/'
#    formatted_ham = get_formatted_hamiltonian(filename, filepath)
#    write_to_textfile(formatted_ham, filename, writepathname)

# Ansatz circuit formatting

# filepath = '../VQEBAB/operator_data/'
# writepathname = '../VQEBAB/operator_data/'
# filename = 'ansatz_symbolic.txt'
# formatted_circuit = get_formatted_circuit(filename, filepath)
# write_to_textfile(formatted_circuit, filename, writepathname)
# formatted_circuit
