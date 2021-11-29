function [condl_u, condl_d, P_s, cond_zero_d, cond_zero_u] = vprobe(TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
    Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd, Tlu_p, Tld_p, Tru_p, Trd_p, lenT, ee, Ln, de, muL, muR, bb_in, voltage)
%VPROBE Summary of this function goes here
%   Detailed explanation goes here

EE = ee;
Fp_lu=zeros(1,Ln*2); Fp_ru=zeros(1,Ln*2); Fp_p=zeros(Ln*2);
Fp_ld=zeros(1,Ln*2); Fp_rd=zeros(1,Ln*2);

Flu_p=zeros(1,Ln*2); Fru_p=zeros(1,Ln*2);
Fld_p=zeros(1,Ln*2); Frd_p=zeros(1,Ln*2);

MM=zeros(2*Ln);

condl_u = zeros(1,lenT);
condl_d = zeros(1,lenT);
P_s = zeros(1,lenT);

for tt=1:1:lenT %loop over temperature
    
    bb = bb_in(tt);
    
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
    
%     condl_u = zeros(nn,lenT);
%     condl_d = zeros(nn,lenT);
%     P_s = zeros(nn,lenT);
    
    %calculating conductance
    
    condl_u(tt) = Il_u/(voltage); %condr_u = Ir_u/(voltage);
    condl_d(tt) = Il_d/(voltage); %condr_d = Ir_d/(voltage);
    
    %         GDL(nn) = sum(TRp_ld(:,1:2*Ln)*(muL-MU(1:2*Ln)))+TRLd_Rd*(muL-muR)+TRLd_Ru*(muL-muR);
    %         GUL(nn) = sum(TRp_lu(:,1:2*Ln)*(muL-MU(1:2*Ln)))+TRLu_Rd*(muL-muR)+TRLu_Ru*(muL-muR);
    
    P_s(tt) = (condl_u(tt)-condl_d(tt))/(condl_u(tt)+condl_d(tt));
    
    %===========Comparison to Guo calculations at zero T================
end

cond_zero_d = (TLd_Rd' + TLd_Ru')*(muL)-(TRd_Ld' + TRu_Ld')*(muR);
cond_zero_u = (TLu_Rd' + TLu_Ru')*(muL)-(TRu_Lu' + TRd_Lu')*(muR);
for ii=1:Ln*2
    cond_zero_d = cond_zero_d + (Tld_p(:,ii))*(muL)-Tp_ld(:,ii)*MU(ii);
    cond_zero_u = cond_zero_u + (Tlu_p(:,ii))*(muL)-Tp_lu(:,ii)*MU(ii);
end

end

