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
% Specify initial tank pressure
% - Calculate initial tank fluid properties
% Find mdot through system for given tank fluid properties
% - Calculate pressure drop through valve, tubing, and regulator for input
%   mass flow
% - Calculate adiabatic temp change of gas due to pressure drop
% - Recalculate rate of discharge with new P & T at orifice
% - Iterate on mdot until input and orifice mass flows match
% Calculate new tank properties after time step
% - Specify time step
% - Calculate fluid mass discharged
% - Calculate new density in tank
% - Calculate P & T from density, assuming isentropic expansion
% - Repeat

%% Initial tank parameters
V = 0.049;% m^3, Standard K cylinder volume is 49 L
A = pi*0.002^2;% Choked orifice area for 1100 psi upstream, 0.1081 kg/s
Ptank0 = convpres(2000,"psi","Pa");% Pa
Ttank0 = 298.15;% K
rhotank0 = PREoS(Methane,"rho",Ptank0,Ttank0);% kg/m^3

%% Initialize system as set of parts, flow coefficients, P and T values
systable = table('Size',[2 6],...
    'VariableTypes',{'string','double','double','double','double','double'},...
    'VariableNames',{'PartName','Cv','P1','P2','T1','T2'});
systable.PartName(1) = "HVMF"; systable.Cv(1) = 0.69;
systable.PartName(2) = "RGMF"; systable.Cv(2) = 0.3;
systable.P2(end) = Ptank0; systable.T2(end) = Ttank0;% Initial orifice conditions

%% Iterate as tank pressure drops
Ptank = Ptank0;
Ttank = Ttank0;% TODO: Add time-dependence & tank draining
rhotank = rhotank0;
tstep = 0.1;% sec
while true

% Iterate to find dP through all parts, and mdot at orifice for given tank P & T
% Create objective function based on feedsystemmdot as function only of
% mdot to iteratively find actual system mdot
system = @(mdot) feedsystemmdot(systable, Methane, Ptank, Ttank, A, mdot);
% Find a bracket of mdot values for which the flow rate delta crosses zero
[mdotll, mdotul] = bracket_sign_change(system,0,0.1);
% Find the "root" of the system via bisection to find mass flow input that
% corresponds to orifice mass flow.
[mdotll, mdotul] = bisection(system,mdotll,mdotul);
mdot_opt = mean([mdotll,mdotul]);
% Get updated system table for 'optimal' mass flow rate
[dmdot, mdot, systable] = feedsystemmdot(systable, mdot_opt, Methane, Ptank, Ttank, A);

% Check if blowdown is finished (orifice is unchoked)
if systable.P2(end) < 101325*((Methane.gam+1)/2)^(Methane.gam/(Methane.gam-1))
    break
end

% New Tank Properties
mtank = rhotank*V - mdot*tstep;
rhotank = mtank/V;
% Objective function comparing real-gas pressure to isentropic expansion
% pressure, both as function of temperature. Use objective function to find
% tank pressure after time step.
Pfn = @(T) PREoS(Methane,"P",T,rhotank) - Ptank*(T/Ttank)^(Methane.gam/(Methane.gam-1));
[Tll, Tul] = bracket_sign_change(Pfn,Ttank*0.9,Ttank);
[Tll, Tul] = bisection(Pfn,Tll,Tul);
Ttank = mean([Tll,Tul]);
Ptank = PREoS(Methane,"P",Ttank,rhotank);

end