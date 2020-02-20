% Gas Feed System Simulation

close all; clear all;

%% Load in fluid properties
FluidDatabase

%% Load in Pc values & associated params
PcTrends = readtable('PcTrends1inthroat.csv');
PcTrends = [repelem({0},length(PcTrends{1,:}));PcTrends];
% plot(PcTrends.mdotF,PcTrends.Pc)
% Linear Fit to get Pc from mdot
mdotO2Pc = @(mdot) interp1(PcTrends.mdotO,PcTrends.Pc,mdot);
Pc2mdotF = @(Pc) interp1(PcTrends.Pc,PcTrends.mdotF,Pc);


%% Build system as table of parts, parameters, P & T values

% Blank System Table
make_systab = @(n_parts) table('Size',[n_parts 12],...
    'VariableTypes',[repelem({'string'},2),repelem({'double'},9),{'string'}],...
    'VariableNames',{'PartName','Type','Cv','RegP2','RegDroop','Cd','A','P1','P2','T1','T2','Choked'});


% % MAIN INJECTOR METHANE
% medium = Methane;
% inj_Pc_targ = 10e5;% Pa, 10 bar
% inj_mdot_targ = 0.07945;% kg/s
% 
% RGMF_P2_reg = 50e5;% Pa
% RGMF_droop = 20349856;% Pa/(kg/s), Regulator outlet pressure droop from Tescom 26-2064D24A270 (catalog flow tables)
% orif_D = 3.85;% mm, Flow Control Orifice Size
% orif_A = pi*(orif_D/(2*1000))^2;% m^2
% orif_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% % A_injF = orifice_size(Methane,Pc2mdot(Pc_target),0.4,0.61,41.8e5,224);% m^2
% injF_A = 2.75e-5;% m^2
% injF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% 
% n_inj_pt = 6;
% inj_sys = make_systab(n_inj_pt);
% inj_sys.Choked = repelem("N",n_inj_pt)';
% inj_sys{1,1:2}=["RGMF","regulator"];inj_sys.Cv(1)=0.3;inj_sys.RegP2(1)=RGMF_P2_reg;inj_sys.RegDroop(1)=RGMF_droop;% Tescom 26-2095TA470AN
% inj_sys{2,1:2}=["BVMF1","valve"];inj_sys.Cv(2)=6.0;inj_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
% inj_sys{3,1:2}=["ORMF","orifice"];inj_sys.Cd(3)=orif_Cd;inj_sys.A(3)=orif_A;% Flow Control Orifice
% inj_sys{4,1:2}=["CKMF","valve"];inj_sys.Cv(4)=1.9;% CheckAll U3CSSTF0.500SS check valve
% inj_sys{5,1:2}=["BVMF2","valve"];inj_sys.Cv(5)=6.0;inj_sys.A(5)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
% inj_sys{6,1:2}=["inj","orifice"];inj_sys.Cd(6)=injF_Cd;inj_sys.A(6)=injF_A;% Main Injector Methane Orifices
% 
% 
% % IGNITER METHANE
% medium = Methane;
% ig_Pc_targ = convpres(200,"psi","Pa");% Pa
% ig_mdot_targ = 0.00110279054076074;% kg/s
% 
% RGIF_P2_reg = 32e5;% Pa
% RGIF_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
% meter_Cv = 0.0098;% Metering Valve Cv<=0.03
% igF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% igF_A = 7.59E-07;% m^2, From igniter spreadsheet
% 
% n_ig_pt = 6;
% ig_sys = make_systab(n_ig_pt);
% ig_sys.Choked = repelem("N",n_ig_pt)';
% ig_sys{1,1:2}=["RGIF","regulator"];ig_sys.Cv(1)=0.1147;ig_sys.RegP2(1)=RGIF_P2_reg;ig_sys.RegDroop(1)=RGIF_droop;% Victor SR4J
% ig_sys{2,1:2}=["BVIF","valve"];ig_sys.Cv(2)=6.0;ig_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
% ig_sys{3,1:2}=["NVIF","valve"];ig_sys.Cv(3)=meter_Cv;% Swagelok SS-4MG2-MH Flow Metering Valve, Vernier Handle
% ig_sys{4,1:2}=["CKIF","valve"];ig_sys.Cv(4)=1.9;% CheckAll U3CSSTF.500SS check valve
% ig_sys{5,1:2}=["SVIF","valve"];ig_sys.Cv(5)=0.04;ig_sys.A(5)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
% ig_sys{6,1:2}=["ig","orifice"];ig_sys.Cd(6)=igF_Cd;ig_sys.A(6)=igF_A;% Igniter Methane Orifice


% IGNITER OXYGEN
medium = Oxygen;
ig_Pc_targ = convpres(200,"psi","Pa");% Pa
ig_mdot_targ = 0.00496;% kg/s

RGIF_P2_reg = 46.7e5;% Pa
RGIF_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
meter_Cv = 0.023;% Metering Valve Cv<=0.03
igF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
igF_A = 2.3E-06;% m^2, From igniter spreadsheet

n_ig_pt = 7;
ig_sys = make_systab(n_ig_pt);
ig_sys.Choked = repelem("N",n_ig_pt)';
ig_sys{1,1:2}=["HVIO","valve"];ig_sys.Cv(1)=0.69;% Sherwood GV Cylinder Valve
ig_sys{2,1:2}=["RGIO","regulator"];ig_sys.Cv(2)=0.1147;ig_sys.RegP2(2)=RGIF_P2_reg;ig_sys.RegDroop(2)=RGIF_droop;% Victor SR4J
ig_sys{3,1:2}=["BVIO","valve"];ig_sys.Cv(3)=6.0;ig_sys.A(3)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
ig_sys{4,1:2}=["NVIO","valve"];ig_sys.Cv(4)=meter_Cv;% Swagelok SS-4MG2-MH Flow Metering Valve, Vernier Handle
ig_sys{5,1:2}=["CKIO","valve"];ig_sys.Cv(5)=1.9;% CheckAll U3CSSTF.500SS check valve
ig_sys{6,1:2}=["SVIO","valve"];ig_sys.Cv(6)=0.04;ig_sys.A(6)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
ig_sys{7,1:2}=["ig","orifice"];ig_sys.Cd(7)=igF_Cd;ig_sys.A(7)=igF_A;% Igniter Ox Orifice


systab = ig_sys;
n_parts = n_ig_pt;
% systab = inj_sys;
% n_parts = n_inj_pt;


%% Simulation Setup

% Termination time, time step
t_step = 0.25;% sec
t_stop = 15;% sec
n_steps = t_stop/t_step+1;


% Initial tank parameters
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
Ptank0 = convpres(2400,"psi","Pa");% Pa
Ttank0 = 273.15+25;% K
rhotank0 = PREoS(medium,"rho",Ptank0,Ttank0);% kg/m^3


% Assemble logging table from part names, to logtab P&T at outlet of each
% part. t, mdot, and tank params are explicitly included for storage
varnames = {'t','mdot','Ptank','Ttank'};
for i=1:n_parts
    varnames{2*i+3} = sprintf('%s_P2',systab.PartName(i));
    varnames{2*i+4} = sprintf('%s_T2',systab.PartName(i));
end
vartypes = repelem({'double'},length(varnames));
logtab = table('Size',[n_steps length(varnames)],...
    'VariableTypes',vartypes,'VariableNames',varnames);
logtab.t(1) = 0;


% Initialize tank conditions
Ptank = Ptank0;
Ttank = Ttank0;
rhotank = rhotank0;


% Add waitbar
wbar = waitbar(0,'Running Simulation');


%% Iterate on tank conditions until t_stop
for itr=1:n_steps

% Find mdot for given tank conditions
mdot_step = 1e-2;% kg/s
mdot_coef = 0.5;% Step size refinement coefficient
mdot_tol = 1e-6;% convergence criterion
mdot_test = mdot_step; % Initialize mdot_test
systab.P1(1) = Ptank;
systab.T1(1) = Ttank;
[mdot_sys, systab] = system_mdot(systab, medium, mdot_test, mdot_step, mdot_coef, mdot_tol);
% while true
%     % Propagate P & T through system for mdot_test
%     mdot_out = 1;% Reset mdot_out > mdot_test to avoid choke boolean later
%     for itm=1:n_parts
%         % Get P1, T1 from nearest upstream element
%         if itm==1
%             systab.P1(itm) = Ptank;
%             systab.T1(itm) = Ttank;
%         else
%             systab.P1(itm) = systab.P2(itm-1);
%             systab.T1(itm) = systab.T2(itm-1);
%         end
% 
%         % Calculate P2, T2 from P1, T1, and check for choking
%         % All flow device codes assume isentropic temp relations
%         if systab.Type(itm) == "valve"
%             itm_out = valve(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cv(itm));
%         elseif systab.Type(itm) == "regulator"
%             itm_out = regulator(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cv(itm),systab.RegP2(itm),systab.RegDroop(itm));
%         elseif systab.Type(itm) == "orifice"
%             itm_out = orifice(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cd(itm),systab.A(itm)); 
%         end
% 
%         % If unchoked, add P2, T2 to system table
%         if length(itm_out) > 1
%             systab.P2(itm) = itm_out(1);
%             systab.T2(itm) = itm_out(2);
%         else % if choked, get output mdot & log choke location
%             mdot_out = itm_out;
%             systab.Choked = repelem("N",n_parts)';
%             systab.Choked(itm) = "Y";
%             break
%         end
%     end
% 
%     % Check convergence
%     if mdot_step <= mdot_tol
%         mdot_sys = mdot_test; % Declare system mdot as converged
%         break
%     end
% 
%     % If choked, go back 1 step & refine step size to "sneak up" on choked
%     % condition. If not choked, continue with same step size.
%     if mdot_out < mdot_test
%         mdot_test = mdot_test - mdot_step;
%         mdot_step = mdot_step*0.5;
%     else
%         mdot_test = mdot_test + mdot_step;
%     end
% 
% end


% logtab values from iteration
if itr ~= 1
    logtab.t(itr) = logtab.t(itr-1)+t_step;
end
logtab.mdot(itr) = mdot_sys;
logtab.Ptank(itr) = Ptank;
logtab.Ttank(itr) = Ttank;
for i=1:n_parts
    logtab(itr,varnames(2*i+3)) = systab(i,'P2');
    logtab(itr,varnames(2*i+4)) = systab(i,'T2');
end


% New Tank Properties
mtank = rhotank*V_tank - mdot_sys*t_step;
rhotank = mtank/V_tank;
% Only new tank density is known, so use isentropic expansion w/ EoS to
% find new temp, then use new temp & known density to find P
Pfn = @(T) PREoS(medium,"P",T,rhotank) - Ptank*(T/Ttank)^(medium.gam/(medium.gam-1));
[Tll, Tul] = bracket_sign_change(Pfn,Ttank*0.9,Ttank);
[Tll, Tul] = bisection(Pfn,Tll,Tul);
Ttank = mean([Tll,Tul]);
Ptank = PREoS(medium,"P",Ttank,rhotank);

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
xlabel('Time, s')
ylabel('Igniter Oxygen Mass Flow Rate, kg/s')


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
title('Pressure Drop Through Igniter Oxygen Feed System')
xlabel('Component')
ylabel('Time, s')
zlabel('Pressure, bar')
xticklabels(Plabels)

figure
plot(xs,station_Ps(1,:),'-k','DisplayName','t = 0 s','LineWidth',2)
hold on
plot(xs,interp1(logtab.t,station_Ps,3.33333),'--k','DisplayName','t = 5 s','LineWidth',2)
plot(xs,interp1(logtab.t,station_Ps,6.66667),'-.k','DisplayName','t = 10 s','LineWidth',2)
plot(xs,station_Ps(end,:),':k','DisplayName','t = 15 s','LineWidth',2)
title('Pressure Drop Through Igniter Oxygen Feed System')
xlabel('Component')
ylabel('Pressure, bar')
legend('Location','northeast')
% xticks(xs)
xticklabels(Plabels)
grid on
hold off

% %% Print key outputs
% fprintf('Injector Outlet Pressure: %0.4f bar\n',logtab.ig_P2(end)/1e5)
% fprintf('                    mdot: %0.4f kg/s\n',logtab.mdot(end))
% fprintf('Chamber Pressure at mdot: %0.4f bar\n',mdotO2Pc(logtab.mdot(end))/1e5)
% disp(systab)


%% Print key outputs
fprintf('Igniter mdot: %0.6f kg/s\n',logtab.mdot(1))
fprintf('Igniter Target mdot: %0.6f kg/s\n',ig_mdot_targ)
% % fprintf('Main Combustion Chamber Pressure at mdot: %0.2f bar\n',mdot2Pc(logtab.ig_mdot(1))/1e5)
fprintf('Igniter Outlet Pressure: %0.2f bar\n',logtab.ig_P2(1)/1e5)
fprintf('Igniter Target Pressure: %0.2f bar\n',ig_Pc_targ/1e5)
disp(systab)