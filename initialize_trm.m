function [TLu_Ru, TLu_Rd, TLd_Ru, TLd_Rd,  TRu_Lu, TRu_Ld, TRd_Lu, TRd_Ld, Tp_p, Tp_lu, Tp_ld, Tp_ru, Tp_rd,Tp2_p, Tlu_p, Tld_p, Tru_p, Trd_p] = initialize_trm(Ln, lenE)

%This function initializes the transmission matrices to save memory, some
%of the variable here are not used

TLu_Ru  = sparse(1,lenE);
TLu_Rd  = sparse(1,lenE);
TLd_Ru  = sparse(1,lenE);
TLd_Rd  = sparse(1,lenE);

% DDD NEW 
TRu_Lu  = sparse(1,lenE);
TRu_Ld  = sparse(1,lenE);
TRd_Lu  = sparse(1,lenE);
TRd_Ld  = sparse(1,lenE);

Tp_p = zeros(lenE, 2*Ln, 2*Ln);
Tp_lu = sparse(lenE, 2*Ln);
Tp_ld = sparse(lenE, 2*Ln);
Tp_ru = sparse(lenE, 2*Ln);
Tp_rd = sparse(lenE, 2*Ln);

Tp2_p = zeros(lenE, 2*Ln, 2*Ln);
Tlu_p = sparse(lenE, 2*Ln);
Tld_p = sparse(lenE, 2*Ln);
Tru_p = sparse(lenE, 2*Ln);
Trd_p = sparse(lenE, 2*Ln);


end

