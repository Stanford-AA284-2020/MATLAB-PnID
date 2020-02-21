% Setup file for individual component (neat) property calculations from
% Helmholtz fundamental relations.
% C.F. Edwards, 2/2/08

% Set the fractional tolerance for numerical solutions.
global toler, toler = 1e-6;

% Universal constants.  Note that these may get changed in the data files
% when they are loaded.  (Ru is a good example.)
global Ru, Ru = 8314.34;

% Critical and triple point data.
global Tcrit_i Pcrit_i rcrit_i
global Ttrip_i Ptrip_i rftrip_i rgtrip_i

% Limits of known accuracy
global Tupper_i Tlower_i Pupper_i rupper_i

% Species constants.
global M_i R_i wP_i 

% Ideal gas reference state info.
global Tref_i Pref_i rref_i href_i sref_i

% Fundamental relation data.
global FR_Npoly0 FR_Neinst FR_Nspec0
global FR_N0 FR_t0 FR_gamma0 FR_spec0
global FR_Npoly FR_Nexp FR_Ngaus FR_Nnonan
global FR_N FR_d FR_t FR_c
global FR_eta FR_beta FR_gamma FR_epsilon
global FR_a FR_b FR_NAbeta FR_A FR_B FR_C FR_D

% Permit access to species indices.
global N2 O2 Ar CO2 H2O nH2 CH4

% Set some limit values for allocating storage.  These are not critical;
% they just need to be larger than actual requirements.
N     = 7;      % Max number of species to be loaded
N0    = 19;     % Max number of ideal gas coefficients
NFR   = 56;     % Max number of residual FR coefficients

% Set the index values for each species to be loaded.  This must be done
% before getting the individual species data so that it goes in the right
% place.
N2  = 1 
O2  = 2 
Ar  = 3
CO2 = 4
H2O = 5
nH2 = 6
CH4 = 7
% Make the data storage.
Tcrit_i = zeros(1,N); Pcrit_i = zeros(1,N); rcrit_i = zeros(1,N);
Ttrip_i = zeros(1,N); Ptrip_i = zeros(1,N); 
rftrip_i = zeros(1,N); rgtrip_i = zeros(1,N);
Tupper_i = zeros(1,N); Tlower_i = zeros(1,N); 
Pupper_i = zeros(1,N); rupper_i = zeros(1,N);
M_i = zeros(1,N); R_i = zeros(1,N); wP_i = zeros(1,N);
Tref_i = zeros(1,N); Pref_i = zeros(1,N); rref_i = zeros(1,N); 
href_i = zeros(1,N); sref_i = zeros(1,N);

FR_Npoly0 = zeros(1,N); FR_Neinst = zeros(1,N); FR_Nspec0 = zeros(1,N);
FR_Npoly = zeros(1,N); FR_Nexp = zeros(1,N); 
FR_Ngaus = zeros(1,N); FR_Nnonan = zeros(1,N);
FR_N0 = zeros(N0,N); FR_t0 = zeros(N0,N); 
FR_gamma0 = zeros(N0,N); FR_spec0 = zeros(N0,N);
FR_N = zeros(NFR,N); FR_d = zeros(NFR,N); FR_t = zeros(NFR,N); FR_c = zeros(NFR,N);
FR_eta = zeros(NFR,N); FR_beta = zeros(NFR,N);
FR_gamma = zeros(NFR,N); FR_epsilon = zeros(NFR,N);
FR_a = zeros(NFR,N); FR_b = zeros(NFR,N);
FR_NAbeta = zeros(NFR,N);
FR_A = zeros(NFR,N); FR_B = zeros(NFR,N);
FR_C = zeros(NFR,N); FR_D = zeros(NFR,N);

% Get the properties data.
Get_Tegeler_Ar_Data
Get_Span_N2_Data
Get_Schmidt_O2_Data
Get_Span_CO2_Data
Get_Wagner_H2O_Data
Get_Jacobsen_nH2_Data
Get_Jacobsen_CH4_Data