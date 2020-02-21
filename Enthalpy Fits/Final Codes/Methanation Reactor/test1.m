%% This script solves the coupled differential equations for concentration
close all; clear all; clc;
tspan = [0.02, 2];

global P;
global velocity;
global rho_cat;
global Q;
global Ru; 
%global T;
T = 523; %Reactor temperature
T_init = T;
Q = 0;
P = 3e6;
velocity = 0.0006;
rho_cat = 4500; % kg/m^3
Ru = 8.314;
y0 = [T; 0.01; 0.2; 0; 0.78; 0.01];
% y0 = [T; 0.5; 0; 0; 0.001; 0.5];
%y0 = [T; 0.211541267927211;   0.307660397158147;   0.009581509608962;   0.423203068473957;   0.048013756831723];
y0 = y0*P/(Ru*T); y0(1,1) = T;
[t,y] = ode15s('methanation_mk9', tspan, y0);
len = size(y,1)

T_final = y(len,1);
y_final = y(len,:);
xs = y_final(1,2:6);
xs = xs*(Ru*T_final/P);
xs = xs/sum(xs)


% T = y(:,1);
% CH4_final = y(:,4);
% CH4_final = CH4_final.*(T);
% CH4_final = CH4_final*Ru/P;
% H2O_final = y(:,6).*(T); H2O_final = H2O_final*Ru/P;
% H2_final = y(:,5).*(T); H2_final = H2_final*Ru/P;
% CO2_final = y(:,3).*(T); CO2_final = CO2_final*Ru/P;
% CO_final = y(:,2).*(T); CO_final = CO_final*Ru/P;

T = []; CO_final = []; CO2_final = []; CH4_final = [];
H2_final = []; H2O_final = [];

for i = 1: size(y,1)
    T_temp = y(i,1);
    CO = y(i,2)*T_temp*Ru/P;
    CO2 = y(i,3)*T_temp*Ru/P;
    CH4 = y(i,4)*T_temp*Ru/P;
    H2 = y(i,5)*T_temp*Ru/P;
    H2O = y(i,6)*T_temp*Ru/P;
    n_tot = CO + CO2 + CH4 + H2 + H2O;
    T = [T; T_temp];
    CO_final = [CO_final; CO/n_tot];
    CO2_final = [CO2_final; CO2/n_tot];
    CH4_final = [CH4_final; CH4/n_tot];
    H2_final = [H2_final; H2/n_tot];
    H2O_final = [H2O_final; H2O/n_tot];
end


figure(1);
plot(t,CH4_final, 'r-'); hold on; plot(t,H2O_final, 'b-'); plot(t,H2_final, 'k-'); plot(t,CO_final, 'g-'); plot(t,CO2_final, 'y-');
xlabel('Distance along the catalyst bed (m)'); ylabel('Mole Fractions'); 
title(['Initial T = ', num2str(T_init), 'K, P = ', num2str(P/1e5), 'bar']);
str = {['Velocity - \it' num2str(velocity) ' m/s '],['Catalyst Density - \it' num2str(rho_cat) ' kg/m^3'],[ 'Heat Input - \it' num2str(Q) ' J/kg-cat']};
txtt = text(0.6,0.4,str); txtt.FontSize = 10.5;
legend('Methane', 'Water', 'Hydrogen', 'Carbon Monoxide', 'Carbon Dioxide');
plotfixer;

figure(2);
plot(t,CH4_final, 'r-');
xlabel('Distance along the catalyst bed (m)'); ylabel('Mole Fractions'); 
title(['Initial T = ', num2str(T_init), 'K, P = ', num2str(P/1e5), 'bar']);
plotfixer;

figure(3);
plot(t,T, 'r-');
xlabel('Distance along the catalyst bed (m)'); ylabel('Temperature (K)'); 
title('Temperature profile within the reactor');
plotfixer;

% %%Equilbrium
% gas = IdealGasMix('gasification_small.xml');
% N     = nSpecies(gas);
% iCO  = speciesIndex(gas,'CO');
% iH2O = speciesIndex(gas,'H2O');
% iCO2  = speciesIndex(gas,'CO2');
% iH2   = speciesIndex(gas,'H2');
% iCH4 = speciesIndex(gas,'CH4');
% xgas = zeros(1,N);
% %CO = []; CO2 =[]; CH4 = []; H2O = []; H2 = []; Tf = [];
% %Equilibrium for Methanation
% % T1 = 273; 
% % P1 = 4e8;
% T1 = 273; 
% P1 = 4e8;
% xgas(iCO) = 0.01;
% xgas(iCO2) = 0.2; 
% xgas(iH2) = 0.78;
% xgas(iH2O) = 0.01;
% gas1 = gas;
% set(gas1, 'T', T1, 'P', P1, 'X', xgas);
% 
% 
% gas2 = equilibrate(gas1, 'HP');
% xgas2 = moleFractions(gas2);
% CO_eq = xgas2(iCO);
% CO2_eq = xgas2(iCO2);
% CH4_eq = xgas2(iCH4);
% H2_eq = xgas2(iH2);
% H2O_eq = xgas2(iH2O);
% 
% figure(1);
% hold on;
% plot([0.02 2], [CO_eq CO_eq], '--g');
% plot([0.02 2], [CO2_eq CO2_eq],'--y');
% plot([0.02 2], [CH4_eq CH4_eq], '--r');
% plot([0.02 2], [H2_eq H2_eq],'--k');
% plot([0.02 2], [H2O_eq H2O_eq], '--b');
% text(0.6, 0.2, {' \it Dotted Lines : Equilibrium concentrations '});
% plotfixer;
