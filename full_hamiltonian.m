function [hh_NNplus] = full_hamiltonian(length2)
%FULL_HAMILTONIAN Summary of this function goes here
%   Detailed explanation goes here
    Ln = length2;
    hh_NNplus = zeros(Ln);
    
    Inputs
    
%     for jj=1:1:Ln %make hamiltonian matrix spin up
%         for ii=1:1:Ln
%             if jj~=ii
%                 hh_NNplus(ii, jj)=tunneling_decay(jj);
%             end
%         end
%     end 
    
    for nn=1:1:Ln
        for jj=1:1:Ln-nn
            
            lj = sqrt((2*R*sin(jj*delta_phi/2))^2+(jj*delta_h)^2);
            theta_j = cosh(2*R*sin(jj*delta_phi/2)/lj);
            phi_nj_p = delta_phi*(nn+jj)+delta_phi*nn;
            phi_nj_n = delta_phi*(nn+jj)-delta_phi*nn;
            
            hh_NNplus(nn, nn+jj)=tunneling_decay(nn,lj);
            hh_NNplus(nn+jj, nn)=tunneling_decay(nn,lj);

        end
    end
end

