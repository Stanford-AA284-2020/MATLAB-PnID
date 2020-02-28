% Code for computing state evolution of the cryogenic run tank
close all; clear all; clc;
% Provide access to support files via the Matlab path.
addpath 'Fundamental Relation Files' 
addpath 'Fundamental Relation Data'
addpath 'Setup Files' 
addpath 'Property Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n**************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

load Results.mat;
time = 0:0.05:10;

% LOx props
T_injector_up = downstream_T(:,3);
P_injector_up = downstream_P(:,3);
d_Ox = 0.00127;
A_ox = pi*(d_Ox^2)/4;

% Surface Tension Correlation
A = 0.0378427336225679;
B = -0.000280199822845628;
C = -0.00000009855225489719;
D = 0.00000000194328853049;

% CH4 props
slope_mdot = (0.0803-0.0792)/10;
slope_density = (12.98 - 12.06)/10;
mdot_fuel = 0.0792:slope_mdot*0.05:0.0803;
density_fuel = 12.06:slope_density*0.05:12.98;
A_fuel = 1E-5;
d_CH4 = ((d_Ox^2)+(4*(1E-5)/pi))^0.5

for i = 1:1:length(time);
    density_Ox(1,i) = rl_iTP(O2, T_injector_up(i,1), P_injector_up(i,1));
    V_exit_Ox(1,i) = mdot(1,i)./(10*density_Ox(1,i)*A_ox);
    J_exit_Ox(1,i) = density_Ox(1,i).*(V_exit_Ox(1,i)^2);
    V_exit_fuel(1,i) = mdot_fuel(1,i)./(10*density_fuel(1,i)*A_fuel);
    J_exit_fuel(1,i) = density_fuel(1,i).*(V_exit_fuel(1,i)^2);
    J_ratio(1,i) = J_exit_fuel(1,i)./J_exit_Ox(1,i);
    Core_Length(1,i) = 38.97811 - 9.99161.*J_ratio(1,i);
    Surface_Tension_Ox(1,i) = A + (B.*(T_injector_up(i,1))) + (C.*(T_injector_up(i,1)^2)) + (D.*(T_injector_up(i,1)^3));
    We(1,i) = density_fuel(1,i).*((V_exit_fuel(1,i) - V_exit_Ox(1,i)).^2).*d_Ox./Surface_Tension_Ox(1,i);
    Flame_Spread_Angle(1,i) = 4.03964 + (0.00466*We(1,i)) - ((1.20615E-7)*(We(1,i)^2));
end

figure(1);
clf;
hold on;
plot(time, J_ratio, 'k-');
xlabel('Time [s]'); ylabel('J');
title('Momentum Flux Ratio');
plotfixer;

figure(2);
clf;
hold on;
plot(time,Core_Length.*(d_Ox*100), 'k-');
xlabel('Time [s]'); ylabel('Length [cm]');
title('Intact Core Length');
plotfixer;

figure(3);
clf;
hold on;
plot(time, Flame_Spread_Angle, 'k-');
xlabel('Time [s]'); ylabel('Degree');
title('Jet Spread Angle');
plotfixer;

    
    