de =1e-3; DD=5;  % grid parameters for energy integration
ee=-DD:de:DD; lenE=length(ee);

Temperature = [200:10:300]; %Temperature in kelvin

ee1=1.60217662e-19; % charge of electron in C, conversion factor from Joule (J) to eV
kb=1.3806485e-23; %Boltzmann constant in J/K  
KB=kb/ee1; %Boltzmann constant in eV/K

TT=Temperature.*KB; % temperature in eV (k_B*T)
bb_in=1./TT; %inverse temperature

epn=0; tn=1; Len=[1:15];
GaL=1; % coupling of first site to L
GaR=1;  % coupling energy of last site to R
GaP=0.06*tn; %dephasing strength

voltage=0.05; %voltage bias on electrodes

muL = 0+voltage/2; muR = 0-voltage/2;

% Spin orbit parameters
delta_phi = 5*pi/9; %change in phi angle
phi_0 = 0; % starting phi angle
theta = 0.66; %theta angle
so = 0.12*tn; % spin orbit factor

%Beyond nearest neighbor parameters
tn_1 = tn; %transmission value first site to second site
len_n = [4.1 5.8 5.1 6.2 8.9 10]; %array of distances of n site to the first site, must match N sites
len_c = 0.9; %exponential decay term
len_1 = 4.1 ; %distance from first site to second site


lenN = length(Len);
lenT = length(TT);
