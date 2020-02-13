%% Gas tank blowdown simulation

% Tank Blowdown Method used in:
%   Haque et al., Rapid depressurization of pressure vessels, J. Loss Prev.
%       Process Ind., Vol 3, January 1990
% - Select Pressure Decrement
% - Perform isentropic flash on each zone
% - Calculate rate of discharge through choke
% - Calculate time step duration
% - Calculate fluid mass discharged
% - Calculate heat transfer coefficients for each zone
% - Perform energy and mass balances over contents of each zone and energy
%   balance over vessel wall
% - If depressurization is complete, stop; otherwise repeat

% Method here:
% - Specify initial static pressure

% - Calculate rate of discharge through orifice with static pressure
% - Calculate pressure drop through valve, tubing, and regulator
% - Calculate adiabatic temp change of gas due to pressure drop
% - Recalculate rate of discharge with pressure drop
% - Repeat until pressure drops & orifice mass flow match

% - Specify time step & Calculate fluid mass discharged
% - Perform isentropic expansion of gas in cylinder to find
% - Repeat until P = Pa

gam = 1.32;% Methane, STP
R = 8314.46;% J/kmol-K, Gas Constant
Rs = R/16.043;% J/kg-K, Specific Gas Constant, Methane
V = 0.049;% m^3, Standard K cylinder volume is 49 L
A = pi*0.001635^2;% Choked orifice area for 1100 psi upstream, 0.1081 kg/s
P0 = convpres(2000,"psi","Pa");% Pa
T0 = 298.15;% K
rho0 = PREoS("Methane","rho",P0,T0);% kg/m^3
tol = 0.1;

Pt = P0;
Tt = t0;
while delta > tol
    mdot = gam/((gam+1)/2)^((gam+1)/(2*(gam-1)))*(Pt*A/sqrt(gam*Rs*Tt));
    