tic %Starts timer

Inputs

calctype = 1;  %calctype, 1 = voltage probe, 2 = dephasing prove

condl_u = zeros(lenN,lenT);
condl_d = zeros(lenN,lenT);
P_s = zeros(lenN,lenT);

for nn=1:lenN

    Ln = Len(nn);
    [hh_up, hh_down] = make_hamiltonian(epn, tn, Ln); %makes hamiltonian for spin up and spin down
    [HH_NN] = full_hamiltonian(Ln); %makes spin hamiltonian for tunneling terms beyond nearest neightbour
    
    
    HH = blkdiag(hh_up+HH_NN,hh_down+HH_NN); %joins the spin up and down hamiltonian together
    
    if so ~= 0
        bloc_spin = make_spinHamNN(theta,phi_0,Ln,delta_phi,so); %makes the SOC hamiltonian contributions
        HH = HH + bloc_spin; %joining hamiltonian and the SOC contribution
        disp(HH)
    end

    [Gamma_LU, Gamma_RU, Gamma_LD, Gamma_RD, Gamma_P, gamma] = make_gamma(GaL, GaR, GaP, Ln); %makes all the gamma matrices

    HH = effective_ham(HH, gamma); %adds gamma contributions to the hamiltonian
    
    [TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
    Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd, Tlu_p, Tld_p, Tru_p, Trd_p] = transmission_alt(Ln, lenE, HH, Gamma_P, Gamma_LU, Gamma_LD, Gamma_RU, Gamma_RD, ee);

% =====================================================================voltage probe stuff

    if calctype == 1
        
        [condl_u_new, condl_d_new, P_s_new, cond_zero_d, cond_zero_u] = vprobe(TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
        Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd, Tlu_p, Tld_p, Tru_p, Trd_p, lenT, ee, Ln, de, muL, muR, bb_in, voltage);
        condl_u(nn,:)=condl_u_new;
        condl_d(nn,:)=condl_d_new;
        P_s(nn,:)=P_s_new;
        
    end
    
%     clear Fld_p Fld_rd Fld_ru Flu_p Flu_rd Flu_ru Fp2_p Fp_ld Fp_lu Fp_p Fp_rd Fp_ru

    if nn == lenN
    
    cond_zero_d = cond_zero_d/(muL-muR); %are factors needed here?
    cond_zero_u = cond_zero_u/(muL-muR);
    
    P_zero = (cond_zero_u-cond_zero_d)./(cond_zero_u+cond_zero_d); %spin polarization
    
    end
    
    %==========================================================================DEPHASING PROBE CALCULATIONS
    
    if calctype == 2
        
        %Building the MM matrix
        
        MM= zeros(lenE, Ln*2, Ln*2);
        
        for aa=1:1:Ln*2
            for aap=1:1:Ln*2
                MM(:,aa,aap) = -TRp_p(:,aa,aap);
            end
            MM(1:lenE,aa,aa)=0;
            for aap=1:1:Ln*2
                if aap~=aa
                    MM(:,aa,aa)=TRp_p(:,aa,aap) + MM(:,aa,aa);
                end
            end
            %Is this correct, is the direction of transmission correct?
            MM(:,aa,aa)=MM(:,aa,aa)+TRp2_ld(:,aa)+TRp2_rd(:,aa)+TRp2_lu(:,aa)+TRp2_ru(:,aa);
        end
        
        %building the v matrix
        
        vec=zeros(lenE,Ln*2); fp=zeros(lenE,Ln*2);
        Ilp=zeros(1,Ln*2);    Irp=zeros(1,Ln*2);
        
        EE = ee;
        
        fl=1./(exp(bb.*(EE-muL))+1); fr=1./(exp(bb.*(EE-muR))+1);
        
        %is the setup of vec correct?
        
        for aa=1:1:Ln*2
            vec(:,aa)=TRp2_lu(:,aa).*fl' + TRp2_ru(:,aa).*fr'+TRp2_ld(:,aa).*fl'+TRp2_rd(:,aa).*fr';
        end
        
        %Solving for Mn
        for ii=1:1:lenE
            Mn(1:Ln*2,1:Ln*2)=MM(ii,1:Ln*2,1:Ln*2);
            if (rcond(Mn) < 1e-12)
                fp(ii,1:Ln*2)=pinv(Mn)*(vec(ii,1:Ln*2)');
            else
                fp(ii,1:Ln*2)=Mn\(vec(ii,1:Ln*2)');
            end
        end
        
        %Conductances from left lead spin to right lead spin
        IlLu_Ru=sum(TRLu_Ru.*(fl-fr))*de;
        IlLu_Rd=sum(TRLu_Rd.*(fl-fr))*de;
        IlLd_Ru=sum(TRLd_Ru.*(fl-fr))*de;
        IlLd_Rd=sum(TRLd_Rd.*(fl-fr))*de;
        
        Ilru = IlLu_Ru + IlLu_Rd; %is adding this okay????
        Ilrd = IlLd_Ru + IlLd_Rd;
        
        %Calculating
        for aa=1:1:Ln*2
            Ilpu(aa)= sum(TRp2_lu(:,aa).*(fl'-fp(:,aa)))*de;
            Irpu(aa)= sum(TRp2_ru(:,aa).*(fr'-fp(:,aa)))*de;
            Ilpd(aa)= sum(TRp2_ld(:,aa).*(fl'-fp(:,aa)))*de;
            Irpd(aa)= sum(TRp2_rd(:,aa).*(fr'-fp(:,aa)))*de;
        end
        
        Il_u=Ilru+sum(Ilpu);        Ir_u=-Ilru+sum(Irpu);
        Il_d=Ilrd+sum(Ilpd);        Ir_d=-Ilrd+sum(Irpd);
        GLu(nn) = Il_u/voltage; GRu(nn) = Ir_u/voltage;
        GLd(nn) = Il_d/voltage; GRd(nn) = Ir_d/voltage;
        
    end
    
    %saving the data
    
%     clear TRd_Ld TRd_Lu TRu_Ld TRu_Lu Tp2_p Tp_p Tp_ld Tp_lu Tp_rd Tp_ru...
%           Tld_p Tlu_p Trd_p Tru_p
    
end

if calctype == 1
    save data_volt_10sites2 cond_zero_d cond_zero_u condl_d...
        condl_u P_s Temperature de DD epn tn Len GaL GaR GaP voltage delta_phi theta so...
        len_n len_c len_1 P_zero
elseif calctype == 2
    save data_deph
end

clearvars -except data_volt_30sites cond_zero_d cond_zero_u condl_d condl_u...
          Temperature de DD epn tn Len GaL GaR GaP voltage delta_phi theta so...
          len_n len_c len_1 P_s P_zero

toc
