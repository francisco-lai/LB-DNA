function [blc_ul,blc_ur,blc_dl,blc_dr] = make_spinmat(so,a11,a12,a21,a22,ii, blc_ul,blc_ur,blc_dl,blc_dr)

    %This function creates a matrix containing the spin contributions for
    %the hamiltonian

    blc_ul(ii,ii+1) = 1i*so*a11;
    blc_ul(ii+1, ii) = -1i*so*conj(a11); %conj

    blc_dl(ii,ii+1) = 1i*so*a21;
    blc_dl(ii+1, ii) = -1i*so*conj(a12); %conj

    blc_ur(ii,ii+1) = 1i*so*a12;
    blc_ur(ii+1, ii) = -1i*so*conj(a21); %conj

    blc_dr(ii,ii+1) = 1i*so*a22;
    blc_dr(ii+1, ii) = -1i*so*conj(a22); %conj

end
