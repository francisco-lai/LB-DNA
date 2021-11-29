tic %Starts timer

Inputs

calctype = 1;  %calctype, 1 = voltage probe, 2 = dephasing prove


for nn=1:lenN
    
    Ln = Len(nn);
    [hh_up, hh_down] = make_hamiltonian(epn, tn, Ln); %makes hamiltonian for spin up and spin down
    [HH_NN] = full_hamiltonian(Ln); %makes spin hamiltonian for tunneling terms beyond nearest neightbour
    HH = blkdiag(hh_up+HH_NN,hh_down+HH_NN); %joins the spin up and down hamiltonian together
    
    if so ~= 0
        bloc_spin = make_spinHamNN(theta,phi_0,Ln,delta_phi,so); %makes the SOC hamiltonian contributions
        HH = HH + bloc_spin; %joining hamiltonian and the SOC contribution
    end

    [Gamma_LU, Gamma_RU, Gamma_LD, Gamma_RD, Gamma_P, gamma] = make_gamma(GaL, GaR, GaP, Ln); %makes all the gamma matrices

    HH = effective_ham(HH, gamma); %adds gamma contributions to the hamiltonian

    [TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
    Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd,Tp2_p, Tlu_p, Tld_p, Tru_p, Trd_p] = initialize_trm(Ln, lenE); %initializes the size of each array/matrix
   
    % Calculation of transmissions
    
    for ii=1:lenE
        Gr=inv(ee(ii)*eye(2*Ln)-HH); Ga=Gr'; %Green's function  matrices
        for aa=1:1:Ln*2
            %transmission for probes to probes
            for aap=1:1:Ln*2
                Tp_p(ii,aa,aap)=trace(Gamma_P(:,:,aa)*Ga*Gamma_P(:,:,aap)*Gr);  % transmission from probe p' up to probe p
                Tp2_p(ii,aa,aap)=trace(Gamma_P(:,:,aap)*Ga*Gamma_P(:,:,aa)*Gr);  % transmission from probe p up to probe p'
            end
            
            Tlu_p(ii,aa)=trace(Gamma_LU*Ga*Gamma_P(:,:,aa)*Gr);   % p --> LU
            Tld_p(ii,aa)=trace(Gamma_LD*Ga*Gamma_P(:,:,aa)*Gr);   % p --> LD 
            Tru_p(ii,aa)=trace(Gamma_RU*Ga*Gamma_P(:,:,aa)*Gr);   % p --> RU
            Trd_p(ii,aa)=trace(Gamma_RD*Ga*Gamma_P(:,:,aa)*Gr);   % p --> RD

            Tp_lu(ii,aa)=trace(Gamma_P(:,:,aa)*Ga*Gamma_LU*Gr);   % LU --> p
            Tp_ld(ii,aa)=trace(Gamma_P(:,:,aa)*Ga*Gamma_LD*Gr);   % LD --> p
            Tp_ru(ii,aa)=trace(Gamma_P(:,:,aa)*Ga*Gamma_RU*Gr);   % RU --> p
            Tp_rd(ii,aa)=trace(Gamma_P(:,:,aa)*Ga*Gamma_RD*Gr);   % RD --> p 


        end

         TRu_Lu(ii)=trace(Gamma_RU*Ga*Gamma_LU*Gr);  % LU --> RU
         TRu_Ld(ii)=trace(Gamma_RU*Ga*Gamma_LD*Gr);  %  LD --> RU
         TRd_Lu(ii)=trace(Gamma_RD*Ga*Gamma_LU*Gr);  % LU-->RD
         TRd_Ld(ii)=trace(Gamma_RD*Ga*Gamma_LD*Gr);  % LD -->RD

         TLu_Ru(ii)=trace(Gamma_LU*Ga*Gamma_RU*Gr);  % RU -->LU
         TLu_Rd(ii)=trace(Gamma_LU*Ga*Gamma_RD*Gr);  %  RD-->LU
         TLd_Ru(ii)=trace(Gamma_LD*Ga*Gamma_RU*Gr);  %  RU-->LD
         TLd_Rd(ii)=trace(Gamma_LD*Ga*Gamma_RD*Gr);   % RD-->LD

    end
    
    clear hh_down hh_up HH HH_NN

% =====================================================================voltage probe stuff
    if calctype == 1
        for tt=1:1:lenT %loop over temperature
            
            bb = bb_in(tt);
            EE = ee;
            
            Fp_lu=zeros(1,Ln*2); Fp_ru=zeros(1,Ln*2); Fp_p=zeros(Ln*2);
            Fp_ld=zeros(1,Ln*2); Fp_rd=zeros(1,Ln*2);
            
            Flu_p=zeros(1,Ln*2); Fru_p=zeros(1,Ln*2); Fp2_p=zeros(Ln*2);
            Fld_p=zeros(1,Ln*2); Frd_p=zeros(1,Ln*2);
            
            MM=zeros(2*Ln);
            
            dfdE=-bb./(4*cosh(bb*EE/2).^2); %derivative of the fermi function with respect to energy
            
            Flu_ru=-sum(TLu_Ru.*dfdE)*de;  % elements of matrix MM
            Flu_rd=-sum(TLu_Rd.*dfdE)*de;  % elements of matrix MM
            Fld_ru=-sum(TLd_Ru.*dfdE)*de;  % elements of matrix MM
            Fld_rd=-sum(TLd_Rd.*dfdE)*de;  % elements of matrix MM

            Fru_lu=-sum(TRu_Lu.*dfdE)*de;  % elements of matrix MM
            Frd_lu=-sum(TRd_Lu.*dfdE)*de;  % elements of matrix MM
            Fru_ld=-sum(TRu_Ld.*dfdE)*de;  % elements of matrix MM
            Frd_ld=-sum(TRd_Ld.*dfdE)*de;  % elements of matrix MM
            
            
            for aa=1:Ln*2
                
                Fp_lu(aa)=-sum(Tp_lu(:,aa).*dfdE')*de;
                Fp_ru(aa)=-sum(Tp_ru(:,aa).*dfdE')*de;
                Fp_ld(aa)=-sum(Tp_ld(:,aa).*dfdE')*de;
                Fp_rd(aa)=-sum(Tp_rd(:,aa).*dfdE')*de;

                Flu_p(aa)=-sum(Tlu_p(:,aa).*dfdE')*de;
                Fru_p(aa)=-sum(Tru_p(:,aa).*dfdE')*de;
                Fld_p(aa)=-sum(Tld_p(:,aa).*dfdE')*de;
                Frd_p(aa)=-sum(Trd_p(:,aa).*dfdE')*de;
                
                
                
                for aap=1:1:Ln*2
                    Fp_p(aa,aap)=-sum(Tp_p(:,aa,aap).*dfdE')*de;
                    MM(aa,aap)=-Fp_p(aa,aap);
                end
            end
            
            for aa=1:1:Ln*2
                MM(aa,aa)=0;
                for aap=1:1:Ln*2
                    if aap~=aa
                        MM(aa,aa)=MM(aa,aa)+Fp_p(aap,aa);
                    end
                end
                MM(aa,aa)=MM(aa,aa)+Fp_lu(aa)+Fp_ru(aa)+Fp_ld(aa)+Fp_rd(aa);
            end
            
            vec=zeros(1,Ln*2);
            for aa=1:Ln*2
                vec(aa)=muL*Flu_p(aa)+muR*Fru_p(aa)+muL*Fld_p(aa)+muR*Frd_p(aa); % DDD
            end
            
            MU = MM\(vec');
            
            % DDD NEW
            Il_up=sum(Flu_p(1:Ln*2))*muL;
            Ip_lu=sum(Fp_lu(1:Ln*2)* MU(1:2*Ln));
            Il_u=-(Fru_lu+Frd_lu)*muR  + (Flu_ru+Flu_rd)*muL  +Il_up - Ip_lu;
            
            Il_dp=sum(Fld_p(1:Ln*2))*muL;
            Ip_ld=sum(Fp_ld(1:Ln*2)* MU(1:2*Ln));
            Il_d=-(Fru_ld+Frd_ld)*muR  + (Fld_ru+Fld_rd)*muL  +Il_dp - Ip_ld;
            
            %calculating conductance
            
            condl_u(nn,tt) = Il_u/(voltage); %condr_u = Ir_u/(voltage);
            condl_d(nn,tt) = Il_d/(voltage); %condr_d = Ir_d/(voltage);
            
            %         GDL(nn) = sum(TRp_ld(:,1:2*Ln)*(muL-MU(1:2*Ln)))+TRLd_Rd*(muL-muR)+TRLd_Ru*(muL-muR);
            %         GUL(nn) = sum(TRp_lu(:,1:2*Ln)*(muL-MU(1:2*Ln)))+TRLu_Rd*(muL-muR)+TRLu_Ru*(muL-muR);
            
            P_s(nn,tt) = (condl_u(nn,tt)-condl_d(nn,tt))/(condl_u(nn,tt)+condl_u(nn,tt));
            
            %===========Comparison to Guo calculations at zero T================
            cond_zero_d = (TLd_Rd' + TLd_Ru')*(muL)-(TRd_Ld' + TRu_Ld')*(muR);
            cond_zero_u = (TLu_Rd' + TLu_Ru')*(muL)-(TRu_Lu' + TRd_Lu')*(muR);
            for ii=1:Ln*2
                cond_zero_d = cond_zero_d + (Tld_p(:,ii))*(muL)-Tp_ld(:,ii)*MU(ii);
                cond_zero_u = cond_zero_u + (Tlu_p(:,ii))*(muL)-Tp_lu(:,ii)*MU(ii);
            end
        end
    end
    
    clear Fld_p Fld_rd Fld_ru Flu_p Flu_rd Flu_ru Fp2_p Fp_ld Fp_lu Fp_p Fp_rd Fp_ru
    
    cond_zero_d = cond_zero_d/(muL-muR); %are factors needed here?
    cond_zero_u = cond_zero_u/(muL-muR);
    
    P_zero = (cond_zero_u-cond_zero_d)./(cond_zero_u+cond_zero_d); %spin polarization
    
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
    
    clear TRd_Ld TRd_Lu TRu_Ld TRu_Lu Tp2_p Tp_p Tp_ld Tp_lu Tp_rd Tp_ru...
          Tld_p Tlu_p Trd_p Tru_p
    
    if calctype == 1
        save data_volt_15sites cond_zero_d cond_zero_u condl_d...
             condl_u P_s Temperature de DD epn tn Len GaL GaR GaP voltage delta_phi theta so...
             len_n len_c len_1 P_zero
    elseif calctype == 2 
        save data_deph
    end
    
end

clearvars -except data_volt_30sites cond_zero_d cond_zero_u condl_d condl_u...
          Temperature de DD epn tn Len GaL GaR GaP voltage delta_phi theta so...
          len_n len_c len_1 P_s P_zero

toc
