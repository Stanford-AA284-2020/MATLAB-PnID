function [P2, P3, P4, P5, T2, T3, T4, T5] = Downstream_P(ispecies, Pmax, mdot)
% Return pressures at different stations for a given initial pressure and
% massflow rate

global Tcrit_i rcrit_i Ttrip_i Tupper_i
global toler
persistent inumcrit Tnumcrit rnumcrit Pnumcrit

% Running through the feed sequence

% % Tank Initial Pressure
% Pmax = 1000; % psi
% Pmax = Pmax*6894.76; % psi to Pa conversion

% Injector Plate (Stations 4 - 5)
% P5 = 350; % psi: Chamber Pressure
P5 = 1.5E6; % Pa
T5 = Saturation_iP(ispecies, P5); % K (Can be varied)
% mdot = 0.21; % kg/s
Cd_5 = 0.6; % Assumption
Area_5 = pi*(0.004^2)/4; % Assume diameter to be 1 cm
[P4, T4] = upstream_Valve(ispecies, P5, T5, Cd_5, Area_5, mdot, Pmax);

% Ball Valve (Stations 3 - 4)
Cd_4 = 0.3; % Assumption
Area_4 = pi*(0.008^2)/4; % Assume diameter to be 1 cm
[P3, T3] = upstream_Valve(ispecies, P4,T4,Cd_4,Area_4,mdot,Pmax);

% Check Valve (Stations 2 - 3)
Cd_3 = 0; % Assumption
Cv_3  = 1.9/15850.3;
Area_3 = pi*(0.0088392^2)/4; % Assume diameter to be 1 cm

if Cd_3 == 0;
    Q_3 = mdot/rl_iTP(ispecies, T3, P3);
    dP23 = (Q_3/Cv_3)^2;
    P2 = P3 + dP23;
    h3 = h_irT(ispecies, rl_iTP(ispecies, T3, P3), T3);
    h2 = h3;
    [r2 T2] = rT_ihP(ispecies, h2, P2);
else
    [P2, T2] = upstream_Valve(ispecies, P3,T3,Cd_3,Area_3,mdot,Pmax);
end

