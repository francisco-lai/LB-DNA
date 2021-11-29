function [t_n] = tunneling_decay(jj)
%l_1 is thedistance for site 1 to two

Inputs

if jj > size(len_n,2)
    t_n = 0;
elseif jj <= size(len_n,2)
    l_n = len_n(jj);
    
    t_n = tn_1*exp(-(l_n-len_1)/len_c);
end

end

