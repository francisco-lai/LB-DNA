function [TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
          Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd, Tlu_p, Tld_p, Tru_p, Trd_p] = transmission_alt(Ln, lenE, HH, Gamma_P, Gamma_LU, Gamma_LD, Gamma_RU, Gamma_RD, ee)
      
%TRANSMISSION Summary of this function goes here
%   Detailed explanation goes here

[TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld,...
    Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd,Tp2_p, Tlu_p, Tld_p, Tru_p, Trd_p] = initialize_trm(Ln, lenE); %initializes the size of each array/matrix

newLn= 2*Ln;

parfor ii=1:lenE
    
    [Gr, Ga] = make_greenf(ee, ii, Ln, HH); %Green's function  matrices
    
    GmLUGa = Gamma_LU*Ga;
    GmLDGa = Gamma_LD*Ga;
    GmRUGa = Gamma_RU*Ga;
    GmRDGa = Gamma_RD*Ga;
    GmLUGr = Gamma_LU*Gr;
    GmLDGr = Gamma_LD*Gr;
    GmRUGr = Gamma_RU*Gr;
    GmRDGr = Gamma_RD*Gr;  
%     GmGa = Gamma_P(:,:,:)*Ga;
%     GmGr = Gamma_P(:,:,:)*Gr;
%     GmGr2 = Gamma_P(:,:,:)*Gr;
%     GmGa2 = Gamma_P(:,:,:)*Ga;


    for aa=1:1:newLn
        %transmission for probes to probes
        GmGa = Gamma_P(:,:,aa)*Ga;
        GmGr = Gamma_P(:,:,aa)*Gr;
        for aap=1:1:newLn
            GmGr2 = Gamma_P(:,:,aap)*Gr;
%             GmGa2 = Gamma_P(:,:,aap)*Ga;
            Tp_p(ii,aa,aap)=GmGa(:).'*reshape(GmGr2.',[],1);  % transmission from probe p' up to probe p
%             Tp2_p(ii,aa,aap)=GmGa2(:).'*reshape(GmGr.',[],1);  % transmission from probe p up to probe p'
        end
        Tlu_p(ii,aa)=GmLUGa(:).'*reshape(GmGr.',[],1);   % p --> LU
        Tld_p(ii,aa)=GmLDGa(:).'*reshape(GmGr.',[],1);   % p --> LD
        Tru_p(ii,aa)=GmRUGa(:).'*reshape(GmGr.',[],1);   % p --> RU
        Trd_p(ii,aa)=GmRDGa(:).'*reshape(GmGr.',[],1);   % p --> RD
        
        Tp_lu(ii,aa)=GmGa(:).'*reshape(GmLUGr.',[],1);   % LU --> p
        Tp_ld(ii,aa)=GmGa(:).'*reshape(GmLDGr.',[],1);   % LD --> p
        Tp_ru(ii,aa)=GmGa(:).'*reshape(GmRUGr.',[],1);   % RU --> p
        Tp_rd(ii,aa)=GmGa(:).'*reshape(GmRDGr.',[],1);   % RD --> p
        
        
    end
    
    TRu_Lu(ii)=GmRUGa(:).'*reshape(GmLUGr.',[],1);  % LU --> RU
    TRu_Ld(ii)=GmRUGa(:).'*reshape(GmLDGr.',[],1);  %  LD --> RU
    TRd_Lu(ii)=GmRDGa(:).'*reshape(GmLUGr.',[],1);  % LU-->RD
    TRd_Ld(ii)=GmRDGa(:).'*reshape(GmLDGr.',[],1);  % LD -->RD
    
    TLu_Ru(ii)=GmLUGa(:).'*reshape(GmRUGr.',[],1);  % RU -->LU
    TLu_Rd(ii)=GmLUGa(:).'*reshape(GmRDGr.',[],1);  %  RD-->LU
    TLd_Ru(ii)=GmLDGa(:).'*reshape(GmRUGr.',[],1);  %  RU-->LD
    TLd_Rd(ii)=GmLDGa(:).'*reshape(GmRDGr.',[],1);   % RD-->LD
    
end
end

