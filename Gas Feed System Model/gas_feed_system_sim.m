% Gas Feed System Simulation

close all; clear all;

%% Load in fluid properties
FluidDatabase

%% Load in Pc values & associated params
PcTrends = readtable('PcTrends1inthroat.csv');
PcTrends = [repelem({0},length(PcTrends{1,:}));PcTrends];
% Linear Fit to get Pc from mdot
mdotO2Pc = @(mdot) interp1(PcTrends.mdotO,PcTrends.Pc,mdot);
Pc2mdotF = @(Pc) interp1(PcTrends.Pc,PcTrends.mdotF,Pc);


%% Build system as table of parts, parameters, P & T values

% Blank System Table
make_systab = @(n_parts) table('Size',[n_parts 12],...
    'VariableTypes',[repelem({'string'},2),repelem({'double'},9),{'string'}],...
    'VariableNames',{'PartName','Type','Cv','RegP2','RegDroop','Cd','A','P1','P2','T1','T2','Choked'});


% IGNITER OXYGEN
medium = Oxygen;
ig_Pc_targ = convpres(200,"psi","Pa");% Pa
ig_mdot_targ = 0.00496;% kg/s

RGIF_P2_reg = 38e5;% Pa
RGIF_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
meter_Cv = 0.043;% Metering Valve Cv<=0.03
% igF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% igF_A = 1.08E-06;% m^2, From igniter spreadsheet

n_ig_pt = 6;
ig_sys = make_systab(n_ig_pt);
ig_sys.Choked = repelem("N",n_ig_pt)';
ig_sys{1,1:2}=["HVIO","valve"];ig_sys.Cv(1)=0.69;% Sherwood GV Cylinder Valve
ig_sys{2,1:2}=["RGIO","regulator"];ig_sys.Cv(2)=0.1147;ig_sys.RegP2(2)=RGIF_P2_reg;ig_sys.RegDroop(2)=RGIF_droop;% Victor SR4J
ig_sys{3,1:2}=["SVIO","valve"];ig_sys.Cv(3)=0.04;ig_sys.A(6)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
ig_sys{4,1:2}=["NVIO","valve"];ig_sys.Cv(4)=meter_Cv;% Swagelok SS-4MG2-MH Flow Metering Valve, Vernier Handle
ig_sys{5,1:2}=["CKIO","valve"];ig_sys.Cv(5)=1.9;% CheckAll U3CSSTF.500SS check valve
ig_sys{6,1:2}=["BVIO","valve"];ig_sys.Cv(6)=0.04;ig_sys.A(6)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
% ig_sys{7,1:2}=["ig","orifice"];ig_sys.Cd(7)=igF_Cd;ig_sys.A(7)=igF_A;% Igniter Ox Orifice


% NITROGEN PURGE
% medium = Nitrogen;
% RGN_P2_reg = convpres(1000,"psi","Pa");% Pa
% RGN_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
% meter_Cv = 0.008;% Medium-Flow Metering Valve Cv<=0.03
% 
% n_n2_pt = 5;
% n2_sys = make_systab(n_n2_pt);
% n2_sys.Choked = repelem("N",n_n2_pt)';
% n2_sys{1,1:2}=["HVN","valve"];n2_sys.Cv(1)=0.69;% Sherwood GV Cylinder Valve
% n2_sys{2,1:2}=["RGN","regulator"];n2_sys.Cv(2)=0.1147;n2_sys.RegP2(2)=RGN_P2_reg;n2_sys.RegDroop(2)=RGN_droop;% Victor SR4J
% n2_sys{3,1:2}=["SVN","valve"];n2_sys.Cv(3)=0.04;n2_sys.A(3)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
% n2_sys{4,1:2}=["NVN","valve"];n2_sys.Cv(4)=meter_Cv;% Swagelok SS-4MG2-MH Flow Metering Valve, Vernier Handle
% n2_sys{5,1:2}=["CKNXX","valve"];n2_sys.Cv(5)=1.9*4;% CheckAll U3CSSTF.500SS check valve


systab = ig_sys;
n_parts = n_ig_pt;
% systab = inj_sys;
% n_parts = n_inj_pt;
% systab = n2_sys;
% n_parts = n_n2_pt;


% Initial tank parameters
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
Ptank = convpres(2400,"psi","Pa");% Pa
Ttank = 273.15+25;% K
rhotank = PREoS(medium,"rho",Ptank,Ttank);% kg/m^3


%% Simulation Setup
% Temperature correlation
% temp_interp = "isothermal";% Neglect Temperature Changes
temp_interp = "adiabatic";% Isenthalpic (Adiabatic)
% temp_interp = "isentropic";% Maximum Temperature Change


% Termination time, time step
t_step = 0.25;% sec
t_stop = 10;% sec
n_steps = t_stop/t_step+1;


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


%% Iterate on tank conditions until t_stop
wbar = waitbar(0,'Running Simulation');
for itr=1:n_steps

% Find mdot for given tank conditions
mdot_step = 1e-4;% kg/s
mdot_coef = 0.5;% Step size refinement coefficient
mdot_tol = 1e-8;% convergence criterion
mdot_test = mdot_step; % Initialize mdot_test
systab.P1(1) = Ptank;
systab.T1(1) = Ttank;
% [mdot_sys, systab] = system_mdot(systab, medium, mdot_test, mdot_step, mdot_coef, mdot_tol);
while true
    % Propagate P & T through system for mdot_test
    choked = false;% Reset mdot_out > mdot_test to avoid choke boolean later
    for itm=1:n_parts
        % Get P1, T1 from nearest upstream element
        if itm==1
            systab.P1(itm) = Ptank;
            systab.T1(itm) = Ttank;
        else
            systab.P1(itm) = systab.P2(itm-1);
            systab.T1(itm) = systab.T2(itm-1);
        end

        % Calculate P2, T2 from P1, T1, and check for choking
        % All flow device codes assume isentropic temp relations
        if systab.Type(itm) == "valve"
            itm_out = valve(medium,...
                mdot_test,...
                systab.P1(itm),...
                systab.T1(itm),...
                systab.Cv(itm),...
                temp_interp);
        elseif systab.Type(itm) == "regulator"
            itm_out = regulator(medium,...
                mdot_test,...
                systab.P1(itm),...
                systab.T1(itm),...
                systab.Cv(itm),...
                systab.RegP2(itm),...
                systab.RegDroop(itm),...
                temp_interp);
        elseif systab.Type(itm) == "orifice"
            itm_out = orifice(medium,...
                mdot_test,...
                systab.P1(itm),...
                systab.T1(itm),...
                systab.Cd(itm),...
                systab.A(itm),...
                temp_interp); 
        end

        % If unchoked, add P2, T2 to system table
        if length(itm_out) > 1
            systab.P2(itm) = itm_out(1);
            systab.T2(itm) = itm_out(2);
        else % if choked, get output mdot & log choke location
            choked = true;
            systab.Choked = repelem("N",n_parts)';
            systab.Choked(itm) = "Y";
            break
        end
    end

    % Check convergence
    if mdot_step <= mdot_tol
        mdot_sys = mdot_test; % Declare system mdot as converged
        break
    end

    % If choked, go back 1 step & refine step size to "sneak up" on choked
    % condition. If not choked, continue with same step size.
    if choked
        mdot_step = mdot_step*0.5;
        mdot_test = max([mdot_test - mdot_step 0]);
    else
        mdot_test = mdot_test + mdot_step;
    end

end


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


%% Plot tank conditions, mass flow rates
figure
% subplot(2,4,1)
% sgtitle('Methane Tank Fluid Conditions & Subsystem Mass Flow Rates')
plot(logtab.t,logtab.Ptank/1e5,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Pressure, bar')
title('Oxygen Tank Pressure')

figure
% subplot(2,4,5)
plot(logtab.t,logtab.Ttank,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Temperature, K')
title('Oxygen Tank Temperature')

figure
plot(logtab.t,logtab.mdot,'k','LineWidth',2)
% ylim([min(logtab.mdot)*0.9,max(logtab.mdot)*1.1])
title('Igniter Oxygen Mass Flow Rate')
xlabel('Time, s')
xlim([0 logtab.t(end)])
ylabel('Mass Flow Rate, kg/s')

% figure
% subplot(2,2,1)
% sgtitle('Tank and Orifice Fluid Conditions During Blowdown')
% plot(logtab.t,logtab.Ptank/1e5,'k','LineWidth',2)
% grid on
% xlabel('Time, s')
% ylabel('Tank Pressure, bar')
% 
% subplot(2,2,3)
% plot(logtab.t,logtab.Ttank,'k','LineWidth',2)
% grid on
% xlabel('Time, s')
% ylabel('Tank Temperature, K')

% subplot(2,2,[2,4])% Show mdot large
% plot(logtab.t,logtab.mdot,'k','LineWidth',2)
% grid on
% ylim([logtab.mdot(end)*0.99,logtab.mdot(1)*1.01])
% xlabel('Time, s')
% ylabel('Igniter Oxygen Mass Flow Rate, kg/s')


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
title('Pressure Drop Through Igniter Oxygen Feed')
xlabel('Component')
ylabel('Time, s')
zlabel('Pressure, bar')
xticklabels(Plabels)

figure
plot(xs,station_Ps(1,:),'-k','DisplayName','t = 0 s','LineWidth',2)
hold on
plot(xs,interp1(logtab.t,station_Ps,logtab.t(end)*1/3),'--k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*1/3),'LineWidth',2)
plot(xs,interp1(logtab.t,station_Ps,logtab.t(end)*2/3),'-.k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*2/3),'LineWidth',2)
plot(xs,station_Ps(end,:),':k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)),'LineWidth',2)
title('Pressure Drop Through Igniter Oxygen Feed')
xlabel('Component')
ylabel('Pressure, bar')
legend('Location','northeast')
xticks(xs)
xticklabels(Plabels)
grid on
hold off

% %% Print key outputs
% fprintf('Injector Outlet Pressure: %0.4f bar\n',logtab.ig_P2(end)/1e5)
% fprintf('                    mdot: %0.4f kg/s\n',logtab.mdot(end))
% fprintf('Chamber Pressure at mdot: %0.4f bar\n',mdotO2Pc(logtab.mdot(end))/1e5)
% disp(systab)


%% Print key outputs
fprintf('mdot: %0.6f kg/s\n',logtab.mdot(1))
% fprintf('Igniter Target mdot: %0.6f kg/s\n',ig_mdot_targ)
% % fprintf('Main Combustion Chamber Pressure at mdot: %0.2f bar\n',mdot2Pc(logtab.ig_mdot(1))/1e5)
fprintf('Outlet Pressure: %0.2f bar\n',logtab{1,end-1}/1e5)
% fprintf('Igniter Target Pressure: %0.2f bar\n',ig_Pc_targ/1e5)
disp(systab)