%% Gas Tank Blowdown Simulation

FluidDatabase % Load in fluid properties

% Method used in:
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

%% Termination Conditions, Time Step
mdot_target = 0.1081;% kg/s
t_step = 0.05;% sec
t_stop = 15;% sec
p_choke = 101325*chokeratio(Methane.gam);% Pa

%% Iterate as tank pressure drops
Ptank = Ptank0;
Ttank = Ttank0;% TODO: Add time-dependence & tank draining
rhotank = rhotank0;
n_steps = t_stop*t_step+1;
store = table('Size',[n_steps 8],...
    'VariableTypes',{'double','double','double','double','double','double','double','double'},...
    'VariableNames',{'t','Ptank','HVMF_P2','RGMF_P2','Ttank','HVMF_T2','RGMF_T2','mdot'});
store.t(1) = 0;
itr = 1;
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
[dmdot, mdot, systable] = feedsystemmdot(systable, Methane, Ptank, Ttank, A, mdot_opt);

% Store values from iteration
if itr ~= 1
    store.t(itr) = store.t(itr-1)+t_step;
end
store.Ptank(itr) = Ptank;
store.HVMF_P2(itr) = systable.P2(1);
store.RGMF_P2(itr) = systable.P2(2);
store.Ttank(itr) = Ttank;
store.HVMF_T2(itr) = systable.T2(1);
store.RGMF_T2(itr) = systable.T2(2);
store.mdot(itr) = mdot;
itr = itr + 1;

% Check if blowdown is finished (mdot, time, orifice pressure criteria)
if store.mdot(itr-1)<=mdot_target || store.t(itr-1)>=t_stop || store.RGMF_P2(itr-1)<p_choke
    break
end

% New Tank Properties
mtank = rhotank*V - mdot*t_step;
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

%%
figure
subplot(2,3,1)
sgtitle('Tank and Orifice Fluid Conditions During Blowdown')
plot(store.t,store.Ptank/1e5,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Tank Pressure, bar')

subplot(2,3,2)
plot(store.t,store.RGMF_P2/1e5,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Orifice Upstream Pressure, bar')

subplot(2,3,[3,6])
plot(store.t,store.mdot,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Mass Flow Rate, kg/s')

subplot(2,3,4)
plot(store.t,store.Ttank,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Tank Temperature, K')

subplot(2,3,5)
plot(store.t,store.RGMF_T2,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Orifice Upstream Temperature, K')