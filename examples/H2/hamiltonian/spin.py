from openfermion import FermionOperator

def spin_operator(n_tot):

    '''
    S^2 operator: S^2 * |psi> = s(s+1) * |psi> if |psi> is an eigenstate of the
    electronic Hamiltonian, because [S^2, H_el], at least if there are no
    magnetic fields


    Args: n_tot = number of spatial orbitals
    '''
    spin_op = FermionOperator()
    for p in range(n_tot):
        for q in range(n_tot):
            # Term 1: S+S-
            p_alpha = 2*p
            p_beta = 2*p+1
            q_alpha = 2*q
            q_beta = 2*q+1

            splus_sminus = (
            FermionOperator(((p_alpha, 1), (p_beta, 0), (q_beta, 1), (q_alpha, 0)))
            )

            spin_op += splus_sminus    

            # Term 2: Sz_p * Sz_q 
            sz_p_sz_q = 1/4 * (
            FermionOperator(((p_alpha, 1), (p_alpha, 0), (q_alpha, 1), (q_alpha, 0))) - 
            FermionOperator(((p_alpha, 1), (p_alpha, 0), (q_beta, 1), (q_beta, 0))) -
            FermionOperator(((p_beta, 1), (p_beta, 0), (q_alpha, 1), (q_alpha, 0))) +
            FermionOperator(((p_beta, 1), (p_beta, 0), (q_beta, 1), (q_beta, 0)))
            )

            spin_op += sz_p_sz_q

        # Term 3: -Sz_p
        sz_p = - 1/2 * (FermionOperator(((p_alpha, 1), (p_alpha, 0))) - FermionOperator(((p_beta, 1), (p_beta, 0))))
        spin_op += sz_p

    return spin_op



