function [hh_up, hh_down] = make_hamiltonian(site_energy, coupling, length)
    % this function creates two hamiltonian matrices (one for spin up and one for spin down) using the a specified
    % length, site energy and the tunneling value and outputs two.
    % hamiltonian matrices
    
    % beyond nearest neighbour implementation, need to include spin
    % interactions too? NN hopping integral is given by
    % t_j=t1*exp(-(l_j-l_1)/l_c), where l_c is the decay constant and l_j
    % is the distance to site l_1
    
    % modularity hierarchy. decay_func->values->hamiltonian?
    %                          ->SOC?
    
    epn = site_energy; tn = coupling; Ln = length;
    hh_up = zeros(Ln) ; hh_down = zeros(Ln);
    
    for jj=1:1:Ln-1 %make hamiltonian matrix spin up
        hh_up(jj,jj)=epn; %; hh_up(jj,jj+1)=tn; hh_up(jj+1,jj)=tn;
    end
    for jj=1:1:Ln-1 %make hamiltonian matrix spin down
        hh_down(jj,jj)=epn; %; hh_down(jj,jj+1)=tn; hh_down(jj+1,jj)=tn;
    end
    
    hh_up(end,end)=epn; % hamiltonian spin up
    hh_down(end,end)= epn; % hamiltonian spin down
end