function [TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
          Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd,Tp2_p, Tlu_p, Tld_p, Tru_p, Trd_p] = transmission(Ln, lenE, HH, Gamma_P, Gamma_LU, Gamma_LD, Gamma_RU, Gamma_RD, ee)
      
%TRANSMISSION Summary of this function goes here
%   Detailed explanation goes here

[TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
    Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd,Tp2_p, Tlu_p, Tld_p, Tru_p, Trd_p] = initialize_trm(Ln, lenE); %initializes the size of each array/matrix

for ii=1:lenE
    
    [Gr, Ga] = make_greenf(ee, ii, Ln, HH); %Green's function  matrices
    
    
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
end

