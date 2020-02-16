%% Gas Feed System Simulation

close all; clear all;

FluidDatabase % Load in fluid properties

% Method:
% Specify system properties & tank parameters
% Calculate mdot through system
%   Calculate dP through each component
%   If component is choked, switch to dp calculation up from chamber with
%   same mdot
%   If no components are choked, check chamber pressure against downstream
%   P2, and reduce/increase mdot until within tolerance
% Repeat until mdot into system 
PcTrends = readtable('/Users/JBR132/Documents/_Stanford/AA284B Propulsion System Design Lab/Propellant-Trades/GCH4_LOX_1inThroat/PcTrends1inthroat.csv');


%% Initialize system as set of parts, flow coefficients, P and T values
systab = table('Size',[6 10],...
    'VariableTypes',[repelem({'string'},2),repelem({'double'},8)],...
    'VariableNames',{'PartName','Type','Cv','P2reg','Cd','A','P1','P2','T1','T2'});
systab{1,1:2}=["HVMF","valve"];systab.Cv(1)=0.69;% Sherwood GV cylinder valve
systab{2,1:2}=["RGMF","regulator"];systab.Cv(2)=0.3;systab.P2reg(2)=76e5;% Tescom 26-2095TA470AN
systab{3,1:2}=["BVMF1","valve"];systab.Cv(3)=6.0;systab.A(3)=pi*(0.281*0.0254)^2;% Swagelok SS-44S6
systab{4,1:2}=["ORMF","orifice"];systab.Cd(4)=0.61;systab.A(4)=pi*(3.35/2000)^2;% Flow Control Orifice
systab{5,1:2}=["CKMF","valve"];systab.Cv(5)=1.9;% CheckAll U3CSSTF0.500SS check valve
systab{6,1:2}=["BVMF1","valve"];systab.Cv(6)=6.0;systab.A(6)=pi*(0.281*0.0254)^2;% Swagelok SS-44S6


%% Initial tank parameters
medium = Methane;
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
D_orifice = 3.38;% mm
A_orifice = pi*((D_orifice/1000)/2)^2;% m^2
Ptank0 = convpres(2000,"psi","Pa");% Pa
Ttank0 = 298.15;% K
rhotank0 = PREoS(medium,"rho",Ptank0,Ttank0);% kg/m^3


%% Termination Conditions, Time Step
mdot_target = 0.079;% kg/s
p_choke = 101325*chokeratio(medium.gam)*10;% Pa
t_step = 0.25;% sec
t_stop = 15;% sec
n_steps = t_stop/t_step+1;


%% Initialize tank conditions
Ptank = Ptank0;
Ttank = Ttank0;
rhotank = rhotank0;

%% Initialize orifice conditions as tank conditions since no flow yet
systab.P2(end) = Ptank0; systab.T2(end) = Ttank0;


%% Assemble logging table from part names, to store P&T at outlet of each
% part. t, mdot, and tank params are explicitly included for storage
n_parts = length(systab.PartName);
varnames = {'t','mdot','Ptank','Ttank'};
for i=1:n_parts
    varnames{2*(i+1)+1} = sprintf('%s_P2',systab.PartName(i));
    varnames{2*(i+1)+2} = sprintf('%s_T2',systab.PartName(i));
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
system = @(mdot) feedsystemmdot(systab, medium, Ptank, Ttank, A_orifice, mdot);
% Find a bracket of mdot values for which the flow rate delta crosses zero
[mdotll, mdotul] = bracket_sign_change(system,1e-3,0.1);
% Find the "root" of the system via bisection to find mass flow input that
% corresponds to orifice mass flow.
[mdotll, mdotul] = bisection(system,mdotll,mdotul);
mdot_opt = mean([mdotll,mdotul]);
% Get updated system table for 'optimal' mass flow rate
outs = feedsystemmdot(systab, medium, Ptank, Ttank, A_orifice, mdot_opt,'all outputs');
dmdot = outs{1};
mdot = outs{2};
systab = outs{3};


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
    store(itr,varnames(2*(i+1)+1)) = systab(i,'P2');
    store(itr,varnames(2*(i+1)+2)) = systab(i,'T2');
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
grid on
xlabel('Time, s')
ylabel('Tank Pressure, bar')

subplot(2,3,2)
plot(store.t,store.RGMF_P2/1e5,'k','LineWidth',2)
grid on
% hold on
% stoptime = interp1(store.RGMF_P2,store.t,p_choke);
% plot(xlim(gca),[p_choke/1e5,p_choke/1e5],':k')
% plot([stoptime,stoptime],ylim(gca),':k')
% hold off
xlabel('Time, s')
ylabel('Orifice Upstream Pressure, bar')

subplot(2,3,[3,6])% Show mdot large
plot(store.t,store.mdot,'k','LineWidth',2)
grid on
ylim([store.mdot(end)*0.99,store.mdot(1)*1.01])
hold on
plot(xlim(gca),[mdot_target,mdot_target],':k')
stoptime = interp1(store.mdot,store.t,mdot_target);
plot([stoptime,stoptime],ylim(gca),':k')
hold off
xlabel('Time, s')
ylabel('Mass Flow Rate, kg/s')

subplot(2,3,4)
plot(store.t,store.Ttank,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Tank Temperature, K')

subplot(2,3,5)
plot(store.t,store.RGMF_T2,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Orifice Upstream Temperature, K')

%% Print mean time derivative of mass flow rate
dmdotdt = mean(rmmissing(diff(store.mdot))/t_step);
fprintf('d/dt mdot = %0.6f kg/s^2\nstop time = %0.2f s\n',dmdotdt,stoptime)