%% Gas Tank Blowdown Simulation

close all; clear all;

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


%% Initialize system as set of parts, flow coefficients, P and T values
systable = table('Size',[2 6],...
    'VariableTypes',{'string','double','double','double','double','double'},...
    'VariableNames',{'PartName','Cv','P1','P2','T1','T2'});
systable.PartName(1) = "HVMF"; systable.Cv(1) = 0.69;
systable.PartName(2) = "RGMF"; systable.Cv(2) = 0.3;


%% Initial tank parameters
medium = Methane;
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
A_orifice = pi*0.002^2;% Choked orifice area for 1100 psi upstream, 0.1081 kg/s
Ptank0 = convpres(2000,"psi","Pa");% Pa
Ttank0 = 298.15;% K
rhotank0 = PREoS(medium,"rho",Ptank0,Ttank0);% kg/m^3


%% Termination Conditions, Time Step
mdot_target = 0.1081;% kg/s
p_choke = 101325*chokeratio(medium.gam);% Pa
t_step = 0.1;% sec
t_stop = 25;% sec
n_steps = t_stop/t_step+1;


%% Initialize tank conditions
Ptank = Ptank0;
Ttank = Ttank0;
rhotank = rhotank0;

%% Initialize orifice conditions as tank conditions since no flow yet
systable.P2(end) = Ptank0; systable.T2(end) = Ttank0;


%% Assemble logging table from part names, to store P&T at outlet of each
% part. t, mdot, and tank params are explicitly included for storage
n_parts = length(systable.PartName);
varnames = {'t','mdot','Ptank','Ttank'};
for i=1:n_parts
    varnames{2*(i+1)+1} = sprintf('%s_P2',systable.PartName(i));
    varnames{2*(i+1)+2} = sprintf('%s_T2',systable.PartName(i));
end
vartypes = cell(1,length(varnames));
vartypes(:) = {'double'};
store = table('Size',[n_steps length(varnames)],...
    'VariableTypes',vartypes,'VariableNames',varnames);
%     'VariableTypes',{'double','double','double','double','double','double','double','double'},...
%     'VariableNames',{'t','mdot','Ptank','Ttank','HVMF_P2','RGMF_P2','HVMF_T2','RGMF_T2'});


%% Initialize time and storage index
store.t(1) = 0;
itr = 1;

%% Iterate as tank pressure drops
while true

    
% Iterate to find dP through all parts, and mdot at orifice for given tank P & T
% Create objective function based on feedsystemmdot as function only of
% mdot to iteratively find actual system mdot
system = @(mdot) feedsystemmdot(systable, medium, Ptank, Ttank, A_orifice, mdot);
% Find a bracket of mdot values for which the flow rate delta crosses zero
[mdotll, mdotul] = bracket_sign_change(system,1e-3,0.1);
% Find the "root" of the system via bisection to find mass flow input that
% corresponds to orifice mass flow.
[mdotll, mdotul] = bisection(system,mdotll,mdotul);
mdot_opt = mean([mdotll,mdotul]);
% Get updated system table for 'optimal' mass flow rate
outs = feedsystemmdot(systable, medium, Ptank, Ttank, A_orifice, mdot_opt,'all outputs');
dmdot = outs{1};
mdot = outs{2};
systable = outs{3};


% Store values from iteration
if itr ~= 1
    store.t(itr) = store.t(itr-1)+t_step;
end
store.mdot(itr) = mdot;
store.Ptank(itr) = Ptank;
store.Ttank(itr) = Ttank;
% store.HVMF_P2(itr) = systable.P2(1);
% store.RGMF_P2(itr) = systable.P2(2);
% store.HVMF_T2(itr) = systable.T2(1);
% store.RGMF_T2(itr) = systable.T2(2);
for i=1:n_parts
    store(itr,varnames(2*(i+1)+1)) = systable(i,'P2');
    store(itr,varnames(2*(i+1)+2)) = systable(i,'T2');
end
itr = itr + 1;


% Check if blowdown is finished
if store.t(itr-1)>=t_stop
    break
end


% New Tank Properties
mtank = rhotank*V_tank - mdot*t_step;
rhotank = mtank/V_tank;
% Objective function comparing real-gas pressure to isentropic expansion
% pressure, both as function of temperature. Use objective function to find
% tank pressure after time step.
Pfn = @(T) PREoS(medium,"P",T,rhotank) - Ptank*(T/Ttank)^(medium.gam/(medium.gam-1));
[Tll, Tul] = bracket_sign_change(Pfn,Ttank*0.9,Ttank);
[Tll, Tul] = bisection(Pfn,Tll,Tul);
Ttank = mean([Tll,Tul]);
Ptank = PREoS(medium,"P",Ttank,rhotank);
end

%% Clear empty table rows if terminated early
for i=1:n_steps
    if store{i,:} == 0
        store{i,:} = NaN;
    end
end
store = rmmissing(store);

%% Plot Simulation Output
figure
subplot(2,3,1)
sgtitle('Tank and Orifice Fluid Conditions During Blowdown')
plot(store.t,store.Ptank/1e5,'k','LineWidth',2)
xlabel('Time, s')
ylabel('Tank Pressure, bar')

subplot(2,3,2)
plot(store.t,store.RGMF_P2/1e5,'k','LineWidth',2)
hold on
plot(xlim(gca),[p_choke,p_choke],':k')
% plot([],ylim(gca),':k')
hold off
xlabel('Time, s')
ylabel('Orifice Upstream Pressure, bar')

subplot(2,3,[3,6])% Show mdot large
plot(store.t,store.mdot,'k','LineWidth',2)
hold on
plot(xlim(gca),[mdot_target,mdot_target],':k')
% plot([],ylim(gca),':k')
hold off
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

%% Print mean time derivative of mass flow rate
dmdotdt = mean(rmmissing(diff(store.mdot))/t_step);
fprintf('d/dt mdot = %0.6f kg/s^2\n',dmdotdt)