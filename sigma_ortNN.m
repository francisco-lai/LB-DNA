function [a11,a12,a21,a22] = sigma_ortNN(phi, theta, ii, jj)
    % This function performs the pauli spin matrix function to obtain the
    % SOC terms used. 
    
    phi_plus = (phi*(jj)+phi*ii)/2;
    
    sig_x = [0 1; 1 0];
    sig_y = [0 -1j; 1j 0];
    sig_z = [1 0; 0 -1];
    sig_orth = (sig_x*sin(phi_plus)*sin(theta)-sig_y*cos(phi_plus)*sin(theta)+sig_z*cos(theta));
    a11 = sig_orth(1,1);
    a12 = sig_orth(1,2);
    a21 = sig_orth(2,1);
    a22 = sig_orth(2,2);
end
