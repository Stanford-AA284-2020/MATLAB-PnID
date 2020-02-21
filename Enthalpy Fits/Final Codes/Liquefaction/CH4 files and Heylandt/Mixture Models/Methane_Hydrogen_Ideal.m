% Mixture model for nH2/CH4/Ar of Lemmon et al.
% J. Phys. Chem. Ref. Data, Vol. 29, p. 331, 2000

Ru = 8314.510;  % J/kmol_K

% Beta is one of the mole fraction exponents for temperature reducing.
Beta(nH2,CH4) = 1;
Beta(CH4,nH2) = Beta(nH2,CH4);
% Beta(nH2,Ar) = 1;
% Beta(Ar,nH2) = Beta(nH2,Ar);
% Beta(CH4,Ar) = 1;
% Beta(Ar,CH4) = Beta(CH4,Ar);

% Phi is the other mole fraction exponents for temperature reducing.
Phi(nH2,CH4) = 1;
Phi(CH4,nH2) = Phi(nH2,CH4);
% Phi(nH2,Ar) = 1;
% Phi(Ar,nH2) = Phi(nH2,Ar);
% Phi(CH4,Ar) = 1;
% Phi(Ar,CH4) = Phi(CH4,Ar);

% Zeta is the binary interaction coefficient for temperature reducing.
Zeta(nH2,CH4) = -0.856350;
Zeta(CH4,nH2) = Zeta(nH2,CH4);
% Zeta(nH2,Ar) = -1.237713;
% Zeta(Ar,nH2) = Zeta(nH2,Ar);
% Zeta(CH4,Ar) = -2.115126;
% Zeta(Ar,CH4) = Zeta(CH4,Ar);

% Xi is the binary interaction coefficient for density reducing.
Xi(nH2,CH4) = -0.00041847;
Xi(CH4,nH2) = Xi(nH2,CH4);
% Xi(nH2,Ar) = -0.00076031;
% Xi(Ar,nH2) = Xi(nH2,Ar);
% Xi(CH4,Ar) =  0.00041232;
% Xi(Ar,CH4) = Xi(CH4,Ar);

% Fmix is the binary interaction coefficient for the excess Helmholtz
% energy.
Fmix(nH2,CH4) = 1;
Fmix(CH4,nH2) = Fmix(nH2,CH4);
% Fmix(nH2,Ar) = 1.121527;
% Fmix(Ar,nH2) = Fmix(nH2,Ar);
% Fmix(CH4,Ar) = 0.597203;
% Fmix(Ar,CH4) = Fmix(CH4,Ar);

% Mixture excess fitting parameters:
% Mix_i = [2 2];
% Mix_j = [-1.4 1.5];
% Mix_N = [-0.00195245 0.00871334];
