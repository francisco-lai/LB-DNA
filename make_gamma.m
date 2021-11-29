function [Gamma_LU, Gamma_RU, Gamma_LD, Gamma_RD, Gamma_P, gamma] = make_gamma(GaL, GaR, GaP, Ln)
%This function takes in the input gamma values for the probes and leads and
%the chain length to build the appropriate matrices used in the calculation
%of the transmissions

% Gamma_L=zeros(Ln);   Gamma_R=zeros(Ln); Gamma_P=zeros(2*Ln,2*Ln,2*Ln);
% block_zero = zeros(Ln);
% Gamma_L(1,1)=GaL; 
% Gamma_R(Ln,Ln)=GaR;
% 
% Gamma_LU = blkdiag(Gamma_L, block_zero); %Gamma left lead spin up
% Gamma_RU = blkdiag(Gamma_R, block_zero); %Gamma right lead spin up
% Gamma_LD = blkdiag(block_zero, Gamma_L); %Gamma left lead spin down
% Gamma_RD = blkdiag(block_zero, Gamma_R); %Gamma right lead spin down
% 
% for aa=1:1:2*Ln
%     Gamma_P(aa,aa,aa)=GaP;
% end
% 
% gamma_leads = Gamma_LU + Gamma_RU + Gamma_LD + Gamma_RD;
% 
% gamma = eye(2*Ln)*GaP+gamma_leads;

Gamma_LU=sparse(Ln*2,Ln*2);   Gamma_RU=sparse(Ln*2,Ln*2); Gamma_P=zeros(2*Ln,2*Ln,2*Ln);
Gamma_LD=sparse(Ln*2,Ln*2);   Gamma_RD=sparse(Ln*2,Ln*2);

Gamma_LU(1,1)= GaL;
Gamma_LD(Ln+1,Ln+1)= GaL;
Gamma_RU(Ln,Ln)= GaR;
Gamma_RD(2*Ln,2*Ln)= GaR;

for aa=1:1:2*Ln
    Gamma_P(aa,aa,aa)=GaP;
end

gamma_leads = Gamma_LU + Gamma_RU + Gamma_LD + Gamma_RD;

gamma = speye(2*Ln)*GaP+gamma_leads;

clear gamma_leads 

end

