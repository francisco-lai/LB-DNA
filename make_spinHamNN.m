function bloc_spin = make_spinHamNN(theta,phi_0,Ln,delta_phi,so)
    % This function takes in the helix geometrical parameters and
    % calculates the spin orbit contributions to the hamiltonian, this
    % function outputs a matrix containing only the spin orbit contribtions
    % to the hamiltonian (i.e the off diagonal elements)

blc_ul = zeros(Ln); blc_ur = zeros(Ln); blc_dl = zeros(Ln); blc_dr = zeros(Ln);
phi = phi_0;

R=2.5;

Inputs

for ii=1:1:Ln
    for jj=ii:1:Ln
        if ii~=jj
            %is the changing theta necessary?
            if jj-ii <= size(len_n,2) && jj <= size(len_n,2)
                theta = acos(2*R*sin((jj-ii)*delta_phi/2)/len_n(jj));
            end
            %disp(theta)
            [a11,a12,a21,a22] = sigma_ortNN(phi,theta,ii,jj); %this function applies the pauli matrice calculations, FLL Looks like theta needs to be variable
            phi = phi_0 + (delta_phi*ii+delta_phi*jj)/2; %updates the new phi used in the new site on the helix
            phi_min = phi_0 + (delta_phi*ii-delta_phi*jj)/2;
            so_n = SOC_decay(so,(jj-ii));
            if ii ~= Ln
                [blc_ul,blc_ur,blc_dl,blc_dr] = make_spinmatNN(so_n, ii, jj, a11,a12,a21,a22,blc_ul,blc_ur,blc_dl,blc_dr,phi_min); %takes the soc values and makes a matrix
            end
        end
    end
end
temp_spin = cat(2,blc_ul,blc_ur);
temp2_spin = cat(2,blc_dl,blc_dr);
bloc_spin = cat(1,temp_spin,temp2_spin); %matrix containing the SOC terms for all the sites in the hamiltonian
end

