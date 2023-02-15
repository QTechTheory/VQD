#!/usr/bin/env python3

from openfermionpyscf import run_pyscf
from pytket.utils import QubitPauliOperator
from openfermion import MolecularData, jordan_wigner, FermionOperator

import pubchempy

#codes from Daniel
from format_for_QuestLink import format_hamiltonian 
from spin import spin_operator 


__doc__="""
Generate Hamiltonian of a molecule, spin, and occupation number operator.

Keys:

Get hamiltonian from compound name

    >ham_from_name("water"):

Generate hamiltonian of He dimer using Jordan-Wigner transformation in questlink format:

    >ham = Hamiltonian([('He',(0,0,0)),('He',(0,0,0.52))])

Store the hamiltonian in easy questlink format

    >ham.jw2questlink("filename.txt")

Save the geometry of the molecule

    >ham.geo2file("filename")

Get spin hamiltonian

    >Spin(n_orbitals)

Get occupation hamiltonian

    >Occupation(n_orbitals)



Generate spin operator using Jordan-Wigner transformation 

"""
class Hamiltonian:
    """
    Generates structure of a hamiltonian.
    """

    def __init__(self, geometry, **kwargs):
        """
        Instantiation of the Hamiltonian object.

        params
            :geometry: list(tuples), list of the atoms and its positions in Angstrom, e.g., [('He',0,0,0),('He',1,0,0)] 
            :kwargs: dist, additional arguments to generate MolecularData. See MolecularData
        """
        #rudimentary check the input
        assert(type(geometry) == list), "Geometry must be a list of tuples"
        #initialization
        self.geometry = geometry
        self.basis = 'sto-3g'
        self.multiplicity = 1
        self.charge = 0
        self.molecule = None
        self.jw_op = None
        self._kwargs = kwargs 
        self.gs_energy = None
        

    def set_molecule(self):
        """
        Set the molecule object as a fermionic operator
        """
        self.molecule = MolecularData(self.geometry, self.basis, self.multiplicity, self.charge, **self._kwargs)
        self.molecule = run_pyscf(self.molecule, run_scf=1,run_fci=0)


    def set_groundstate_energy(self):
        """
        Calculate the groundstate energy and set it to attribute gs_energy
        """
        self.molecule = MolecularData(self.geometry, self.basis, self.multiplicity, self.charge, **self._kwargs)
        self.molecule = run_pyscf(self.molecule, run_scf=1,run_fci=1)
        self.gs_energy = self.molecule.fci_energy


    def set_jw_op(self, compress = True):
        """
        Set attribute op_jw, i.e., the operator in the form of Jordan Wigner

        :compress: (True), compress for smaller number of qubits
        """
        #reset the molecule
        self.set_molecule()

        self.jw_op = jordan_wigner(self.molecule.get_molecular_hamiltonian())
        if compress : 
            self.jw_op.compress()


    def jw2questlink(self, outfile=None):
        """
        Get the Jordan-wigner transformation in the questlink-easy formatting

        :outfile: str, if given, it stores to the corresponding text file.
        """
        if self.jw_op == None :
            self.set_jw_op()
        return format2questlink(self.jw_op, outfile=outfile)


    def geo2file(self, outfile):
        """
        Store geometry of the compount to the file

        :outfile: path of output file
        """
        with open(outfile, 'w') as f: 
            for t in self.geometry:
                f.write(' '.join(str(s) for s in t) + '\n')



class Spin:
    """
    Return S^2 operator: 
    S^2 * |psi> = s(s+1) * |psi> if |psi> is an eigenstate of the
    electronic Hamiltonian, because [S^2, H_el], at least if there are no
    magnetic fields


    params 
    :norbital: Number of spatial orbitals. Usually the number of electrons in the atom. 

    """
    def __init__(self, norbital):
        """
        Give an instance of spin fermion operator FermionOperator format

        params 
        :norbital: Number of spatial orbitals. Usually the number of electrons in the atom. 
        """
        assert type(norbital)==int, "spatial orbital must be integer"
        self.norbital = norbital
        self.fermionop = spin_operator(norbital)  
        self.jw_op = None

    
    def set_jw_op(self, compress = True):
        """
        Set attribute op_jw, i.e., the operator in the form of Jordan Wigner

        :compress: (True), compress for smaller number of qubits
        """
        self.jw_op = jordan_wigner(self.fermionop)
        if compress : 
            self.jw_op.compress()


    def jw2questlink(self,outfile=None):
        """
        Format the spin operator in the questlink format

        outfile = str, if given, it stores to the corresponding text file.
        """
        if self.jw_op == None :
            self.set_jw_op()
        return format2questlink(self.jw_op, outfile=outfile)



class Occupation:
    """
    Return S^2 operator: 
    S^2 * |psi> = s(s+1) * |psi> if |psi> is an eigenstate of the
    electronic Hamiltonian, because [S^2, H_el], at least if there are no
    magnetic fields


    params 
    :norbital: Number of spatial orbitals. Usually the number of electrons in the atom. 

    """
    def __init__(self, norbital):
        """
        Give an instance of spin fermion operator FermionOperator format

        params 
        :norbital: Number of spatial orbitals. Usually the number of electrons in the atom. 
        """
        assert type(norbital)==int, "spatial orbital must be integer"
        self.norbital = norbital
        self.fermionop = self.number_operator() 
        self.jw_op = None


    def number_operator(self):
        '''
        Operator that counts the number of particles
        Args: n_tot = number of spatial orbitals
        '''
        number_op = FermionOperator()
        for i in range(self.norbital * 2): # n_qubits = n_tot * 2
            n_i = FermionOperator(((i, 1), (i, 0)))
            number_op += n_i
        return number_op
    

    def set_jw_op(self, compress = True):
        """
        Set attribute op_jw, i.e., the operator in the form of Jordan Wigner

        :compress: (True), compress for smaller number of qubits
        """
        self.jw_op = jordan_wigner(self.fermionop)
        if compress : 
            self.jw_op.compress()


    def jw2questlink(self, outfile):
        """
        Format the spin operator in the questlink format

        outfile = str, if given, it stores to the corresponding text file.
        """
        if self.jw_op == None :
            self.set_jw_op()
        return format2questlink(self.jw_op, outfile=outfile)




def format2questlink(operator, outfile=None):
    """
    Return the questlink formatted

    params
    :operator: The operator output format from openfermion
    :outfile: str, if given, it stores to the corresponding text file.
    """
    opstr = str(operator).replace('\n','').split('+') 
    parsed = []
    for line in opstr:
        split_line = line.split(' ', 1)
        scalar = split_line[0]
        gates = split_line[1].split(']')[0] + ']'
        parsed += [(scalar, gates)]

    formatted = format_hamiltonian(parsed)
    if outfile != None : 
        assert type(outfile)==str, "outfile is a string path"
        with open(outfile, 'w') as off :
            off.write(formatted)

    return formatted


def ham_from_compound(comp):
    """
    Return hamiltonian objects given compound object from pubchempy 

    :comp: object from pubchempy
    """
    geo=[]
    for at in comp.atoms :
        atd = at.to_dict()
        if comp.coordinate_type == '2d':
            geo+=[(atd['element'],(atd['x'],atd['y'],0))]
        else:
            geo+=[(atd['element'],(atd['x'],atd['y'],atd['z']))]
    ham = Hamiltonian(geo)
    ham.charge = comp.charge

    return ham



def ham_from_name(name):
    """
    Return the hamiltonians object given a chemial name.

    :name: string, chemical name 
    """
    comps=pubchempy.get_compounds(name,'name')
    hams=[]
    for comp in comps:
        hams += [ham_from_compound(comp)]

    return hams


