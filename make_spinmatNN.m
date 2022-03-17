function [blc_ul,blc_ur,blc_dl,blc_dr] = make_spinmatNN(so,ii, jj, a11,a12,a21,a22, blc_ul,blc_ur,blc_dl,blc_dr,phi_min)

    %This function creates a matrix containing the spin contributions for
    %the hamiltonian
    
    phi_minus = phi_min;
    
    blc_ul(ii,ii) = 0;
    blc_ul(ii,ii+jj) = 2i*so*cos(phi_minus)*a11;
    blc_ul(ii+jj,ii) = -2i*so*cos(phi_minus)*conj(a11); %conj
    
    blc_dl(ii,ii) = 0;
    blc_dl(ii,ii+jj) = 2i*so*cos(phi_minus)*a21;
    blc_dl(ii+jj,ii) = -2i*so*cos(phi_minus)*conj(a12); %conj

    blc_ur(ii,ii) = 0;
    blc_ur(ii,ii+jj) = 2i*so*cos(phi_minus)*a12;
    blc_ur(ii+jj,ii) = -2i*so*cos(phi_minus)*conj(a21); %conj

    blc_dr(ii,ii) = 0;
    blc_dr(ii,ii+jj) = 2i*so*cos(phi_minus)*a22;
    blc_dr(ii+jj,ii) = -2i*so*cos(phi_minus)*conj(a22); %conj

end

