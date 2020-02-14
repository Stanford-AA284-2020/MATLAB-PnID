%% Gas Tank Blowdown Simulation

FluidDatabase % Load in fluid properties

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

%% Initial tank parameters
% V = 0.049;% m^3, Standard K cylinder volume is 49 L
% A = pi*0.001635^2;% Choked orifice area for 1100 psi upstream, 0.1081 kg/s
Ptank0 = convpres(2000,"psi","Pa");% Pa
Ttank0 = 298.15;% K
% rho0 = PREoS(Methane,"rho",Ptank0,Ttank0);% kg/m^3

%% Initialize system as set of parts, flow coefficients, P and T values
store = table('Size',[2 6],...
    'VariableTypes',{'string','double','double','double','double','double'},...
    'VariableNames',{'PartName','Cv','P1','P2','T1','T2'});
store.PartName(1) = "HVMF"; store.Cv(1) = 0.69;
store.PartName(2) = "RGMF"; store.Cv(2) = 0.3;
store.P2(end) = Ptank0; store.T2(end) = Ttank0;% Initial orifice conditions

%% Iterate to find dP through all parts, and mdot at orifice
mdotll = 0; mdotul = 0.5;

system = @(mdot) systemmdot(store, mdot, Methane, Ptank, Ttank, A);

% mdot_opt = golden_section_search(system, mdotll, mdotul, 100)

Ptank = Ptank0;
Ttank = Ttank0;

mds = 0.05:0.0001:0.25;
dmds = zeros(length(mds),1);
for i=1:length(mds)
    dmds(i) = systemmdot(store,mds(i),Methane,Ptank,Ttank,A);
end
[dmdopt,idx] = min(abs(dmds));
mdot_opt = mds(idx)
[dmdot, mdot2, store] = systemmdot(store, mdot_opt, Methane, Ptank, Ttank, A)

function [delta_mdot, mdot2, store] = systemmdot(store,mdot,medium,Ptank,Ttank,A)
    for i=1:length(store.PartName)
        if i==1
            store.P1(i) = Ptank;
            store.T1(i) = Ttank;
        else
            store.P1(i) = store.P2(i-1);
            store.T1(i) = store.T2(i-1);
        end
        [store.P2(i), store.T2(i)] = valve(medium,mdot,store.P1(i),store.T1(i),store.Cv(i));
    end
    mdot2 = chokedorifice(A,medium.gam,medium.Mw,store.P2(end),store.T2(end));
    delta_mdot = mdot - mdot2;
end