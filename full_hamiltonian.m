function [hh_NNplus] = full_hamiltonian(length)
%FULL_HAMILTONIAN Summary of this function goes here
%   Detailed explanation goes here
    Ln = length;
    hh_NNplus = zeros(Ln);
    
%     for jj=1:1:Ln %make hamiltonian matrix spin up
%         for ii=1:1:Ln
%             if jj~=ii
%                 hh_NNplus(ii, jj)=tunneling_decay(jj);
%             end
%         end
%     end 
    
    for jj=1:1:Ln
        for ii=1:1:Ln
            if jj~=ii
                hh_NNplus(ii, jj)=tunneling_decay(abs(ii-jj));
            end
        end
    end
end

