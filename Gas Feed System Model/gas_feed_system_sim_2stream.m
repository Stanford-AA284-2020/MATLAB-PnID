% Gas Feed System Simulation for 2 Subsystems (Methane Stream)

close all; clear all;

%% Load in fluid properties
FluidDatabase

%% Load in Pc values & associated params
PcTrends = readtable('PcTrends1inthroat.csv');
PcTrends = [repelem({0},length(PcTrends{1,:}));PcTrends];
% plot(PcTrends.mdotF,PcTrends.Pc)
% Linear Fit to get Pc from mdot
mdot2Pc = @(mdot) interp1(PcTrends.mdotF,PcTrends.Pc,mdot);
Pc2mdot = @(Pc) interp1(PcTrends.Pc,PcTrends.mdotF,Pc);


%% Build system as table of parts, parameters, P & T values

% Blank System Table
make_systab = @(n_parts) table('Size',[n_parts 12],...
    'VariableTypes',[repelem({'string'},2),repelem({'double'},9),{'string'}],...
    'VariableNames',{'PartName','Type','Cv','RegP2','RegDroop','Cd','A','P1','P2','T1','T2','Choked'});

% TANK VALVE
HVF_Cv=0.69;% Sherwood GV cylinder valve
tank_sys = make_systab(1);
tank_sys.Choked = "N";
tank_sys{1,1:2}=["HVF","valve"];tank_sys.Cv=HVF_Cv;


% MAIN INJECTOR
inj_Pc_targ = 10e5;% Pa, 10 bar
inj_mdot_targ = Pc2mdot(inj_Pc_targ);% kg/s

RGMF_P2_reg = 95e5;% Pa
RGMF_droop = 20349856;% Pa/(kg/s), Regulator outlet pressure droop from Tescom 26-2064D24A270 (catalog flow tables)
orif_D = 2.65;% mm, Flow Control Orifice Size
orif_A = pi*(orif_D/(2*1000))^2;% m^2
orif_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% A_injF = orifice_size(Methane,Pc2mdot(Pc_target),0.4,0.61,41.8e5,224);% m^2
injF_A = 0.97e-5;% m^2
injF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit

n_inj_pt = 6;
inj_sys = make_systab(n_inj_pt);
inj_sys.Choked = repelem("N",n_inj_pt)';
inj_sys{1,1:2}=["RGMF","regulator"];inj_sys.Cv(1)=0.3;inj_sys.RegP2(1)=RGMF_P2_reg;inj_sys.RegDroop(1)=RGMF_droop;% Tescom 26-2095TA470AN
inj_sys{2,1:2}=["BVMF1","valve"];inj_sys.Cv(2)=6.0;inj_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
inj_sys{3,1:2}=["ORMF","orifice"];inj_sys.Cd(3)=orif_Cd;inj_sys.A(3)=orif_A;% Flow Control Orifice
inj_sys{4,1:2}=["CKMF","valve"];inj_sys.Cv(4)=1.9;% CheckAll U3CSSTF0.500SS check valve
inj_sys{5,1:2}=["BVMF2","valve"];inj_sys.Cv(5)=6.0;inj_sys.A(5)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
inj_sys{6,1:2}=["inj","orifice"];inj_sys.Cd(6)=injF_Cd;inj_sys.A(6)=injF_A;% Main Injector Methane Orifices


% IGNITER
ig_Pc_targ = convpres(200,"psi","Pa");% Pa
ig_mdot_targ = 0.00110279054076074;% kg/s

RGIF_P2_reg = 95e5;% Pa
RGIF_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
meter_Cv = 0.16;% Metering Valve Cv<=0.16
igF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
igF_A = 7.59E-07;% m^2, From igniter spreadsheet

n_ig_pt = 6;
ig_sys = make_systab(n_ig_pt);
ig_sys.Choked = repelem("N",n_ig_pt)';
ig_sys{1,1:2}=["RGIF","regulator"];ig_sys.Cv(1)=0.1147;ig_sys.RegP2(1)=RGIF_P2_reg;ig_sys.RegDroop(1)=RGIF_droop;% Victor SR4J
ig_sys{2,1:2}=["BVIF","valve"];ig_sys.Cv(2)=6.0;ig_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
ig_sys{3,1:2}=["NVIF","valve"];ig_sys.Cv(3)=meter_Cv;% Swagelok SS-4L2-MH Flow Metering Valve, Vernier Handle
ig_sys{4,1:2}=["CKIF","valve"];ig_sys.Cv(4)=1.9;% CheckAll U3CSSTF.500SS check valve
ig_sys{5,1:2}=["SVIF","valve"];ig_sys.Cv(5)=0.04;ig_sys.A(5)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
ig_sys{6,1:2}=["ig","orifice"];ig_sys.Cd(6)=igF_Cd;ig_sys.A(6)=igF_A;% Igniter Methane Orifice


%% Simulation Setup

% Termination time, time step
t_step = 0.25;% sec
t_stop = 15;% sec
n_steps = t_stop/t_step+1;


% Initial tank parameters
medium = Methane;
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
Ptank = convpres(2000,"psi","Pa");% Pa
Ttank = 273.15+25;% K
rhotank = PREoS(medium,"rho",Ptank,Ttank);% kg/m^3


% Assemble logging table from part names, to logtab P&T at outlet of each
% part. t, mdot, and tank params are explicitly included for storage
varnames = {'t','mdot','inj_mdot','ig_mdot','Ptank','Ttank','HVF_P2','HVF_T2'};
for i=1:n_inj_pt
    varnames{2*i+7} = sprintf('%s_P2',inj_sys.PartName(i));
    varnames{2*i+8} = sprintf('%s_T2',inj_sys.PartName(i));
end
for i=1:n_ig_pt
    varnames{2*i+7+2*n_inj_pt} = sprintf('%s_P2',ig_sys.PartName(i));
    varnames{2*i+8+2*n_inj_pt} = sprintf('%s_T2',ig_sys.PartName(i));
end
vartypes = repelem({'double'},length(varnames));
logtab = table('Size',[n_steps length(varnames)],...
    'VariableTypes',vartypes,'VariableNames',varnames);
logtab.t(1) = 0;


% Make waitbar
wbar = waitbar(0,'Running Simulation');


%% Iterate on tank conditions until t_stop
for itr=1:n_steps

% Find injector mdot for given tank outlet conditions
inj_mdot_step = 1e-2;% kg/s
inj_mdot_coef = 0.5;% Step size refinement coefficient
inj_mdot_tol = 1e-6;% convergence criterion
inj_mdot_test = inj_mdot_step; % Initialize mdot_test

% Find igniter mdot for given tank outlet conditions
ig_mdot_step = 1e-4;% kg/s
ig_mdot_coef = 0.5;% Step size refinement coefficient
ig_mdot_tol = 1e-8;% convergence criterion
ig_mdot_test = ig_mdot_step; % Initialize mdot_test

% Find system mdot by interatively increasing until choked
while true
% Propagate P & T thru tank valve
% tank_mdot_choked = 1;% Reset mdot_choked > mdot_test to avoid choke boolean if not choked
tank_mdot_test = inj_mdot_test+ig_mdot_test;
tank_sys.P1 = Ptank;
tank_sys.T1 = Ttank;
itm_out = valve(medium,tank_mdot_test,tank_sys.P1,tank_sys.T1,tank_sys.Cv);
if length(itm_out) > 1 % If unchoked, store P2, T2
    HVF_P2 = itm_out(1);
    HVF_T2 = itm_out(2);
else
    % If choked, go back 1 step & refine step size to "sneak up" on choked
    % condition. If not choked, continue evaluation to subsystems.
%     tank_mdot_choked = itm_out;
    tank_sys.Choked = "Y";
    inj_mdot_test = inj_mdot_test - inj_mdot_step;
    inj_mdot_step = inj_mdot_step*0.5;
    ig_mdot_test = ig_mdot_test - ig_mdot_step;
    ig_mdot_step = ig_mdot_step*0.5;
    continue
end
    
    
% If tank valve is not choked, propagate P & T thru injector & igniter streams
% MAIN INJECTOR
inj_mdot_choked = 1;% Reset mdot_choked > mdot_test to avoid choke boolean if not choked
for itm=1:n_inj_pt
    % Get P1, T1 from nearest upstream element
    if itm~=1
        inj_sys.P1(itm) = inj_sys.P2(itm-1);
        inj_sys.T1(itm) = inj_sys.T2(itm-1);
    end

    % Calculate P2, T2 from P1, T1, and check for choking
    % All flow device codes assume isentropic temp relations
    if inj_sys.Type(itm) == "valve"
        itm_out = valve(medium,inj_mdot_test,inj_sys.P1(itm),inj_sys.T1(itm),inj_sys.Cv(itm));
    elseif inj_sys.Type(itm) == "regulator"
        itm_out = regulator(medium,inj_mdot_test,inj_sys.P1(itm),inj_sys.T1(itm),inj_sys.Cv(itm),inj_sys.RegP2(itm),inj_sys.RegDroop(itm));
    elseif inj_sys.Type(itm) == "orifice"
        itm_out = orifice(medium,inj_mdot_test,inj_sys.P1(itm),inj_sys.T1(itm),inj_sys.Cd(itm),inj_sys.A(itm)); 
    end

    % If unchoked, add P2, T2 to system table
    if length(itm_out) > 1
        inj_sys.P2(itm) = itm_out(1);
        inj_sys.T2(itm) = itm_out(2);
    else % if choked, get output mdot & log choke location
        inj_mdot_choked = itm_out;
        inj_sys.Choked = repelem("N",n_inj_pt)';
        inj_sys.Choked(itm) = "Y";
        break
    end
end
    
    
% IGNITER
ig_mdot_choked = 1;% Reset mdot_choked > mdot_test to avoid choke boolean if not choked
for itm=1:n_ig_pt
    % Get P1, T1 from nearest upstream element
    if itm~=1
        ig_sys.P1(itm) = ig_sys.P2(itm-1);
        ig_sys.T1(itm) = ig_sys.T2(itm-1);
    end

    % Calculate P2, T2 from P1, T1, and check for choking
    % All flow device codes assume isentropic temp relations
    if ig_sys.Type(itm) == "valve"
        itm_out = valve(medium,ig_mdot_test,ig_sys.P1(itm),ig_sys.T1(itm),ig_sys.Cv(itm));
    elseif ig_sys.Type(itm) == "regulator"
        itm_out = regulator(medium,ig_mdot_test,ig_sys.P1(itm),ig_sys.T1(itm),ig_sys.Cv(itm),ig_sys.RegP2(itm),ig_sys.RegDroop(itm));
    elseif ig_sys.Type(itm) == "orifice"
        itm_out = orifice(medium,ig_mdot_test,ig_sys.P1(itm),ig_sys.T1(itm),ig_sys.Cd(itm),ig_sys.A(itm)); 
    end

    % If unchoked, add P2, T2 to system table
    if length(itm_out) > 1
        ig_sys.P2(itm) = itm_out(1);
        ig_sys.T2(itm) = itm_out(2);
    else % if choked, get output mdot & log choke location
        ig_mdot_choked = itm_out;
        ig_sys.Choked = repelem("N",n_ig_pt)';
        ig_sys.Choked(itm) = "Y";
        break
    end
end

    %% Check convergence
    if inj_mdot_step <= inj_mdot_tol && ig_mdot_step <= ig_mdot_tol
        sys_mdot = tank_mdot_test; % Declare subsystem mdot as converged
        break
    end
    
    if inj_mdot_choked < inj_mdot_test
        % If choked, go back 1 step & refine step size to "sneak up" on choked
        % condition. If not choked, continue with same step size.
        inj_mdot_test = inj_mdot_test - inj_mdot_step;
        inj_mdot_step = inj_mdot_step*0.5;
    else
        inj_mdot_test = inj_mdot_test + inj_mdot_step;
    end
    
    if ig_mdot_choked < ig_mdot_test
        % If choked, go back 1 step & refine step size to "sneak up" on choked
        % condition. If not choked, continue with same step size.
        ig_mdot_test = ig_mdot_test - ig_mdot_step;
        ig_mdot_step = ig_mdot_step*0.5;
    else
        ig_mdot_test = ig_mdot_test + ig_mdot_step;
    end


end


% % logtab values from iteration
% if itr ~= 1
%     logtab.t(itr) = logtab.t(itr-1)+t_step;
% end
% logtab.mdot(itr) = mdot_sys;
% logtab.Ptank(itr) = Ptank;
% logtab.Ttank(itr) = Ttank;
% for i=1:n_parts
%     logtab(itr,varnames(2*i+3)) = systab(i,'P2');
%     logtab(itr,varnames(2*i+4)) = systab(i,'T2');
% end
% 
% 
% % New Tank Properties
% mtank = rhotank*V_tank - mdot_sys*t_step;
% rhotank = mtank/V_tank;
% % Only new tank density is known, so use isentropic expansion w/ EoS to
% % find new temp, then use new temp & known density to find P
% Pfn = @(T) PREoS(medium,"P",T,rhotank) - Ptank*(T/Ttank)^(medium.gam/(medium.gam-1));
% [Tll, Tul] = bracket_sign_change(Pfn,Ttank*0.9,Ttank);
% [Tll, Tul] = bisection(Pfn,Tll,Tul);
% Ttank = mean([Tll,Tul]);
% Ptank = PREoS(medium,"P",Ttank,rhotank);

% update waitbar
    waitbar(itr/n_steps,wbar)
end
close(wbar)


%% Plot Simulation Output as mdot, in & out P & T
figure
subplot(2,2,1)
sgtitle('Tank and Orifice Fluid Conditions During Blowdown')
plot(logtab.t,logtab.Ptank/1e5,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Tank Pressure, bar')

subplot(2,2,3)
plot(logtab.t,logtab.Ttank,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Tank Temperature, K')

subplot(2,2,[2,4])% Show mdot large
plot(logtab.t,logtab.mdot,'k','LineWidth',2)
grid on
ylim([logtab.mdot(end)*0.99,logtab.mdot(1)*1.01])
% hold on
% plot(xlim(gca),[mdot_target,mdot_target],':k')
% stoptime = interp1(logtab.mdot,logtab.t,mdot_target);
% plot([stoptime,stoptime],ylim(gca),':k')
% hold off
xlabel('Time, s')
ylabel('Mass Flow Rate, kg/s')

% subplot(2,3,2)
% plot(logtab.t,logtab.injector_P2/1e5,'k','LineWidth',2)
% grid on
% % hold on
% % stoptime = interp1(logtab.RGMF_P2,logtab.t,p_choke);
% % plot(xlim(gca),[p_choke/1e5,p_choke/1e5],':k')
% % plot([stoptime,stoptime],ylim(gca),':k')
% % hold off
% xlabel('Time, s')
% ylabel('Injector Downstream Pressure, bar')

% subplot(2,3,5)
% plot(logtab.t,logtab.injector_T2,'k','LineWidth',2)
% grid on
% xlabel('Time, s')
% ylabel('Injector Downstream Temperature, K')

%% Plot Simulation Output as station conditions
P_idxs = [3,5:2:length(varnames)];
station_Ps = logtab{:,P_idxs}/1e5;
T_idxs = P_idxs+1;
station_Ts = logtab(:,T_idxs);
xs = 1:length(P_idxs);
ys = logtab.t;
[X, Y] = meshgrid(xs,ys);
labels = logtab.Properties.VariableNames;
Plabels = ["Tank",systab.PartName{:}];

figure
mesh(X,Y,station_Ps)
view(45,15)
title('Pressure Drop Through Main Methane Feed System')
xlabel('Component')
ylabel('Time, s')
zlabel('Pressure, bar')
xticklabels(Plabels)

figure
plot(xs,station_Ps(1,:),'-k','DisplayName','t = 0 s','LineWidth',2)
hold on
plot(xs,station_Ps(find(logtab.t == 5),:),'--k','DisplayName','t = 5 s','LineWidth',2)
plot(xs,station_Ps(find(logtab.t == 10),:),'-.k','DisplayName','t = 10 s','LineWidth',2)
plot(xs,station_Ps(end,:),':k','DisplayName','t = 15 s','LineWidth',2)
title('Pressure Drop Through Main Methane Feed System')
xlabel('Component')
ylabel('Pressure, bar')
legend('Location','northeast')
% xticks(xs)
xticklabels(Plabels)
grid on
hold off

%% Print key outputs
fprintf('Injector Outlet Pressure: %0.4f bar\n',logtab.injector_P2(end)/1e5)
fprintf('                    mdot: %0.4f kg/s\n',logtab.mdot(end))
fprintf('Chamber Pressure at mdot: %0.4f bar\n',mdot2Pc(logtab.mdot(end))/1e5)
disp(systab)