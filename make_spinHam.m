function bloc_spin = make_spinHam(theta,phi_0,Ln,delta_phi,so)
    % This function takes in the helix geometrical parameters and
    % calculates the spin orbit contributions to the hamiltonian, this
    % function outputs a matrix containing only the spin orbit contribtions
    % to the hamiltonian (i.e the off diagonal elements)

blc_ul = zeros(Ln); blc_ur = zeros(Ln); blc_dl = zeros(Ln); blc_dr = zeros(Ln);
phi = phi_0;

for ii=1:1:Ln
    [a11,a12,a21,a22] = sigma_ort(phi,theta,ii); %this function applies the pauli matrice calculations
    phi = phi_0 + delta_phi*ii; %updates the new phi used in the new site on the helix
    if ii ~= Ln
        [aa11,aa12,aa21,aa22] = sigma_ort(phi,theta,ii+1); %calculates the pauli matrix calculations for the next site over (n+1)
        a11 = a11+aa11;
        a12 = a12+aa12;
        a21 = a21+aa21;
        a22 = a22+aa22;
        [blc_ul,blc_ur,blc_dl,blc_dr] = make_spinmat(so, a11,a12,a21,a22,ii,blc_ul,blc_ur,blc_dl,blc_dr); %takes the soc values and makes a matrix
    end
end
temp_spin = cat(2,blc_ul,blc_ur);
temp2_spin = cat(2,blc_dl,blc_dr);
bloc_spin = cat(1,temp_spin,temp2_spin); %matrix containing the SOC terms for all the sites in the hamiltonian
end

