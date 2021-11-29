function [so_n] = SOC_decay(so_1,ii)
Inputs
if ii > size(len_n,2)
    so_n = 0;
elseif ii <= size(len_n,2)  
    l_n = len_n(ii);
    
    so_n = so_1*exp(-(l_n-len_1)/len_c);
end

end
