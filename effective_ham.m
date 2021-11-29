
function hh_vert = effective_ham(hamiltonian_mat, gamma)
    %adds gammas to the hamiltonian to obtain the effective hamiltonian
    hh = hamiltonian_mat;
    hh_vert = hh+1i*gamma/2;
end

