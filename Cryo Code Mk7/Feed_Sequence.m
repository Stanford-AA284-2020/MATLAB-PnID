% Running through the feed sequence
clear all; close all; clc;

% Call the Cryo Run Tank Simulator
Cryo_RunTank_Sim_Mk2;

% Tank Initial Pressure
Pmax = 1000; % psi
Pmax = Pmax*6894.76; % psi to Pa conversion

% Injector Plate (Stations 4 - 5)
% P5 = 350; % psi: Chamber Pressure
P5 = 1E6; % Pa
T5 = Saturation_iP(O2, P5); % K (Can be varied)
mdot = 0.21; % kg/s
Cd_5 = 0.9; % Assumption
Area_5 = pi*(0.0022^2)/4; % Assume diameter to be 1 cm
[P4, T4] = upstream_Valve(O2, P5, T5, Cd_5, Area_5, mdot, Pmax);

% Ball Valve (Stations 3 - 4)
Cd_4 = 0.3; % Assumption
Area_4 = pi*(0.008^2)/4; % Assume diameter to be 1 cm
[P3, T3] = upstream_Valve(O2, P4,T4,Cd_4,Area_4,mdot,Pmax);

% Check Valve (Stations 2 - 3)
Cd_3 = 0; % Assumption
Cv_3  = 1.9/15850.3;
Area_3 = pi*(0.0088392^2)/4; % Assume diameter to be 1 cm

if Cd_3 == 0;
    Q_3 = mdot/rl_iTP(O2, T3, P3);
    dP23 = (Q_3/Cv_3)^2;
    P2 = P3 + dP23;
    h3 = h_irT(O2, rl_iTP(O2, T3, P3), T3);
    h2 = h3;
    [r2 T2] = rT_ihP(O2, h2, P2);
else
    [P2, T2] = upstream_Valve(O2, P3,T3,Cd_3,Area_3,mdot,Pmax);
end

% Cavitating Venturi (Stations 1 - 2)
Ratio_Down_Up = P2./P;

figure(3);
clf;
plot(time, Ratio_Down_Up, 'b-');
xlabel('Time [s]'); ylabel('P_2/P_1');
title('Ratio of pressure across venturi');
plotfixer;