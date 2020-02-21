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


%% Initial tank parameters
medium = Methane;
V_tank = 0.049;% m^3, Standard K cylinder volume is 49 L
Ptank = convpres(2000,"psi","Pa");% Pa
Ttank = 273.15+25;% K
rhotank = PREoS(medium,"rho",Ptank,Ttank);% kg/m^3


%% Build system as table of parts, parameters, P & T values

% Blank System Table
make_systab = @(n_parts) table('Size',[n_parts 12],...
    'VariableTypes',[repelem({'string'},2),repelem({'double'},9),{'string'}],...
    'VariableNames',{'PartName','Type','Cv','RegP2','RegDroop','Cd','A','P1','P2','T1','T2','Choked'});

% TANK & HAND VALVE
HVF_Cv=0.69;% Sherwood GV cylinder valve
tank_sys = make_systab(1);
tank_sys.Choked = "N";
tank_sys{1,1:2}=["HVF","valve"];tank_sys.Cv=HVF_Cv;


% MAIN INJECTOR
inj_Pc_targ = 10e5;% Pa, 10 bar
inj_mdot_targ = Pc2mdot(inj_Pc_targ);% kg/s

RGMF_P2_reg = 50e5;% Pa
RGMF_droop = 20349856;% Pa/(kg/s), Regulator outlet pressure droop from Tescom 26-2064D24A270 (catalog flow tables)
orif_D = 3.85;% mm, Flow Control Orifice Size
orif_A = pi*(orif_D/(2*1000))^2;% m^2
orif_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
% A_injF = orifice_size(Methane,Pc2mdot(Pc_target),0.4,0.61,41.8e5,224);% m^2
injF_A = 2.75e-5;% m^2
injF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit

n_inj_pt = 6;
inj_sys = make_systab(n_inj_pt);
inj_sys.Choked = repelem("N",n_inj_pt)';
inj_sys{1,1:2}=["RGMF","regulator"];inj_sys.Cv(1)=0.3;inj_sys.RegP2(1)=RGMF_P2_reg;inj_sys.RegDroop(1)=RGMF_droop;% Tescom 26-2095TA470AN
inj_sys{2,1:2}=["BVMF1","valve"];inj_sys.Cv(2)=6.0;inj_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
inj_sys{3,1:2}=["ORMF","orifice"];inj_sys.Cd(3)=orif_Cd;inj_sys.A(3)=orif_A;% Mass Flow Metering Orifice
inj_sys{4,1:2}=["CKMF","valve"];inj_sys.Cv(4)=1.9;% CheckAll U3CSSTF0.500SS check valve
inj_sys{5,1:2}=["BVMF2","valve"];inj_sys.Cv(5)=6.0;inj_sys.A(5)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
inj_sys{6,1:2}=["inj","orifice"];inj_sys.Cd(6)=injF_Cd;inj_sys.A(6)=injF_A;% Main Injector Methane Orifices


% IGNITER
ig_Pc_targ = convpres(200,"psi","Pa");% Pa
ig_mdot_targ = 0.00110279054076074;% kg/s

RGIF_P2_reg = 32e5;% Pa
RGIF_droop = 111924205;% Pa/(kg/s), Regulator outlet pressure droop from Victor SR4J
meter_Cv = 0.0098;% Metering Valve Cv<=0.03
igF_Cd = 0.61;% Sharp-edged plate orifice in high-Re limit
igF_A = 7.59E-07;% m^2, From igniter spreadsheet

n_ig_pt = 6;
ig_sys = make_systab(n_ig_pt);
ig_sys.Choked = repelem("N",n_ig_pt)';
ig_sys{1,1:2}=["RGIF","regulator"];ig_sys.Cv(1)=0.1147;ig_sys.RegP2(1)=RGIF_P2_reg;ig_sys.RegDroop(1)=RGIF_droop;% Victor SR4J
ig_sys{2,1:2}=["BVIF","valve"];ig_sys.Cv(2)=6.0;ig_sys.A(2)=pi*(0.281/2*0.0254)^2;% Swagelok SS-44S6
ig_sys{3,1:2}=["NVIF","valve"];ig_sys.Cv(3)=meter_Cv;% Swagelok SS-4MG2-MH Flow Metering Valve, Vernier Handle
ig_sys{4,1:2}=["CKIF","valve"];ig_sys.Cv(4)=1.9;% CheckAll U3CSSTF.500SS check valve
ig_sys{5,1:2}=["SVIF","valve"];ig_sys.Cv(5)=0.04;ig_sys.A(5)=pi*(3/64/2*0.0254)^2;% Parker Skinner 71216SN2FU00N0C111C2 Solenoid Valve
ig_sys{6,1:2}=["ig","orifice"];ig_sys.Cd(6)=igF_Cd;ig_sys.A(6)=igF_A;% Igniter Methane Orifice


%% Simulation Setup
% Temperature correlation
% temp_interp = "isothermal";% Neglect Temperature Changes
temp_interp = "adiabatic";% Isenthalpic (Adiabatic)
% temp_interp = "isentropic";% Maximum Temperature Change


% Termination time, time step
t_step = 0.2;% sec
t_stop = 15;% sec
n_steps = t_stop/t_step+1;


% Assemble logging table from part names, to logtab P&T at outlet of each
% part. t, mdot, and tank params are explicitly included for storage
varnames = {'t','sys_mdot','inj_mdot','ig_mdot','Ptank','Ttank','HVF_P2','HVF_T2'};
inj_itr = @(i) 2*i+7;
ig_itr = @(i) 2*i+7+2*n_inj_pt;
for i=1:n_inj_pt
    varnames{inj_itr(i)} = sprintf('%s_P2',inj_sys.PartName(i));
    varnames{inj_itr(i)+1} = sprintf('%s_T2',inj_sys.PartName(i));
end
for i=1:n_ig_pt
    varnames{ig_itr(i)} = sprintf('%s_P2',ig_sys.PartName(i));
    varnames{ig_itr(i)+1} = sprintf('%s_T2',ig_sys.PartName(i));
end
vartypes = repelem({'double'},length(varnames));
logtab = table('Size',[n_steps length(varnames)],...
    'VariableTypes',vartypes,'VariableNames',varnames);
logtab.t(1) = 0;


%% Iterate on tank conditions until t_stop
wbar = waitbar(0,'Running Simulation');
for itr=1:n_steps

% injector settings
inj_mdot_step = 1e-2;% kg/s
inj_step_coef = 0.5;% Step size refinement coefficient
inj_mdot_tol = 1e-6;% convergence criterion
inj_mdot_test = inj_mdot_step; % Initialize mdot_test

% igniter settings
ig_mdot_step = 1e-4;% kg/s
ig_step_coef = 0.5;% Step size refinement coefficient
ig_mdot_tol = 1e-8;% convergence criterion
ig_mdot_test = ig_mdot_step; % Initialize mdot_test

% Find system mdot by interatively increasing until choked
while true
    % Propagate P & T thru tank valve
    tank_mdot_test = inj_mdot_test+ig_mdot_test;
    tank_sys.P1 = Ptank;
    tank_sys.T1 = Ttank;
    itm_out = valve(medium,...
        tank_mdot_test,...
        tank_sys.P1,...
        tank_sys.T1,...
        tank_sys.Cv,...
        temp_interp);
    if length(itm_out) > 1 % If unchoked, store P2, T2
        tank_sys.Choked = "N";
        tank_sys.P2 = itm_out(1);
        inj_sys.P1(1) = tank_sys.P2;
        ig_sys.P1(1) = tank_sys.P2;
        tank_sys.T2 = itm_out(2);
        inj_sys.T1(1) = tank_sys.T2;
        ig_sys.T1(1) = tank_sys.T2;
    else % if choked, log choke location, step back, and refine step size
        tank_sys.Choked = "Y";
        inj_sys.Choked = repelem("N",n_inj_pt)';
        inj_mdot_test = max([inj_mdot_test - inj_mdot_step 0]);
        inj_mdot_step = inj_mdot_step*inj_step_coef;
        ig_sys.Choked = repelem("N",n_ig_pt)';
        ig_mdot_test = max([ig_mdot_test - ig_mdot_step 0]);
        ig_mdot_step = ig_mdot_step*ig_step_coef;
        continue
    end


    %% If tank valve is not choked, propagate P & T thru injector & igniter streams
    % MAIN INJECTOR
    inj_choked = false;
    for itm=1:n_inj_pt
        % Get P1, T1 from nearest upstream element
        if itm~=1
            inj_sys.P1(itm) = inj_sys.P2(itm-1);
            inj_sys.T1(itm) = inj_sys.T2(itm-1);
        end

        % Calculate P2, T2 from P1, T1, and check for choking
        % All flow device codes assume isentropic temp relations
        if inj_sys.Type(itm) == "valve"
            itm_out = valve(medium,...
                inj_mdot_test,...
                inj_sys.P1(itm),...
                inj_sys.T1(itm),...
                inj_sys.Cv(itm),...
                temp_interp);
        elseif inj_sys.Type(itm) == "regulator"
            itm_out = regulator(medium,...
                inj_mdot_test,...
                inj_sys.P1(itm),...
                inj_sys.T1(itm),...
                inj_sys.Cv(itm),...
                inj_sys.RegP2(itm),...
                inj_sys.RegDroop(itm),...
                temp_interp);
        elseif inj_sys.Type(itm) == "orifice"
            itm_out = orifice(medium,...
                inj_mdot_test,...
                inj_sys.P1(itm),...
                inj_sys.T1(itm),...
                inj_sys.Cd(itm),...
                inj_sys.A(itm),...
                temp_interp); 
        end

        % If unchoked, add P2, T2 to system table
        if length(itm_out) > 1
            inj_sys.P2(itm) = itm_out(1);
            inj_sys.T2(itm) = itm_out(2);
        else % if choked, log choke location
            inj_sys.Choked = repelem("N",n_inj_pt)';
            inj_sys.Choked(itm) = "Y";
            inj_choked = true;
            break % break out of subsystem for loop
        end
    end


    % IGNITER
    ig_choked = false;
    for itm=1:n_ig_pt
        % Get P1, T1 from nearest upstream element
        if itm~=1
            ig_sys.P1(itm) = ig_sys.P2(itm-1);
            ig_sys.T1(itm) = ig_sys.T2(itm-1);
        end

        % Calculate P2, T2 from P1, T1, and check for choking
        % All flow device codes assume isentropic temp relations
        if ig_sys.Type(itm) == "valve"
            itm_out = valve(medium,...
                ig_mdot_test,...
                ig_sys.P1(itm),...
                ig_sys.T1(itm),...
                ig_sys.Cv(itm),...
                temp_interp);
        elseif ig_sys.Type(itm) == "regulator"
            itm_out = regulator(medium,...
                ig_mdot_test,...
                ig_sys.P1(itm),...
                ig_sys.T1(itm),...
                ig_sys.Cv(itm),...
                ig_sys.RegP2(itm),...
                ig_sys.RegDroop(itm),...
                temp_interp);
        elseif ig_sys.Type(itm) == "orifice"
            itm_out = orifice(medium,...
                ig_mdot_test,...
                ig_sys.P1(itm),...
                ig_sys.T1(itm),...
                ig_sys.Cd(itm),...
                ig_sys.A(itm),...
                temp_interp); 
        end

        % If unchoked, add P2, T2 to system table
        if length(itm_out) > 1
            ig_sys.P2(itm) = itm_out(1);
            ig_sys.T2(itm) = itm_out(2);
        else % if choked, log choke location
            ig_sys.Choked = repelem("N",n_ig_pt)';
            ig_sys.Choked(itm) = "Y";
            ig_choked = true;
            break % break out of subsystem for loop
        end
    end

    %% Check convergence
    if inj_mdot_step <= inj_mdot_tol && ig_mdot_step <= ig_mdot_tol
        break
    end
    
    %% Update Step Sizes depending on choke condition
    if inj_choked % If choked, step back & refine step size
        inj_mdot_step = inj_mdot_step*inj_step_coef;
        inj_mdot_test = max([inj_mdot_test - inj_mdot_step 0]);% Prevent negative mdot
    else % else step forward
        inj_mdot_test = inj_mdot_test + inj_mdot_step;
    end
    
    if ig_choked % If choked, step back & refine step size
        ig_mdot_step = ig_mdot_step*ig_step_coef;
        ig_mdot_test = max([ig_mdot_test - ig_mdot_step 0]);% Prevent negative mdot
    else % else step forward
        ig_mdot_test = ig_mdot_test + ig_mdot_step;
    end

end


% logtab values from iteration
if itr ~= 1
    logtab.t(itr) = logtab.t(itr-1)+t_step;
end
logtab.sys_mdot(itr) = tank_mdot_test; % Declare mdot as converged
logtab.inj_mdot(itr) = inj_mdot_test;
logtab.ig_mdot(itr) = ig_mdot_test;
logtab.Ptank(itr) = Ptank;
logtab.Ttank(itr) = Ttank;
logtab.HVF_P2(itr) = tank_sys.P2;
logtab.HVF_T2(itr) = tank_sys.T2;
for i=1:n_inj_pt
    logtab(itr,varnames(inj_itr(i))) = inj_sys(i,'P2');
    logtab(itr,varnames(inj_itr(i)+1)) = inj_sys(i,'T2');
end
for i=1:n_ig_pt
    logtab(itr,varnames(ig_itr(i))) = ig_sys(i,'P2');
    logtab(itr,varnames(ig_itr(i)+1)) = ig_sys(i,'T2');
end


% New Tank Properties
rhotank = (rhotank*V_tank - logtab.sys_mdot(itr)*t_step)/V_tank;
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
title('Methane Tank Pressure')

figure
% subplot(2,4,5)
plot(logtab.t,logtab.Ttank,'k','LineWidth',2)
grid on
xlabel('Time, s')
ylabel('Temperature, K')
title('Methane Tank Temperature')

figure
subplot(1,3,1)% Show mdots large
sgtitle('Methane Mass Flow Rates')
plot(logtab.t,logtab.sys_mdot,'k','LineWidth',2)
ylim([min(logtab.sys_mdot)*0.9,max(logtab.sys_mdot)*1.1])
title('Total')
xlabel('Time, s')
xlim([0 logtab.t(end)])
ylabel('Mass Flow Rate, kg/s')

subplot(1,3,2)% Show mdots large
plot(logtab.t,logtab.inj_mdot,'k','LineWidth',2)
ylim([min(logtab.inj_mdot)*0.9,max(logtab.inj_mdot)*1.1])
title('Main Injector')
xlabel('Time, s')
xlim([0 logtab.t(end)])
% ylabel('Mass Flow Rate, kg/s')

subplot(1,3,3)% Show mdots large
plot(logtab.t,logtab.ig_mdot,'k','LineWidth',2)
ylim([min(logtab.ig_mdot)*0.9,max(logtab.ig_mdot)*1.1])
title('Igniter')
xlabel('Time, s')
xlim([0 logtab.t(end)])
% ylabel('Mass Flow Rate, kg/s')


%% Plot Injector station conditions
inj_P_idx = [5,7,inj_itr(1:n_inj_pt)];
inj_Ps = logtab{:,inj_P_idx}/1e5;
T_idxs = inj_P_idx+1;
inj_Ts = logtab(:,T_idxs);
xs = 1:length(inj_P_idx);
ys = logtab.t;
[X, Y] = meshgrid(xs,ys);
% labels = logtab.Properties.VariableNames;
Plabels = ["Tank","HVF",inj_sys.PartName{:}];

% figure
% mesh(X,Y,inj_Ps)
% view(45,15)
% title('Pressure Drop Through Main Injector Methane Feed')
% xlabel('Component')
% ylabel('Time, s')
% zlabel('Pressure, bar')
% xticks(xs)
% xticklabels(Plabels)

figure
plot(xs,inj_Ps(1,:),'-k','DisplayName','t = 0 s','LineWidth',2)
hold on
plot(xs,interp1(logtab.t,inj_Ps,logtab.t(end)*1/3),'--k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*1/3),'LineWidth',2)
plot(xs,interp1(logtab.t,inj_Ps,logtab.t(end)*2/3),'-.k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*2/3),'LineWidth',2)
plot(xs,inj_Ps(end,:),':k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)),'LineWidth',2)
title('Pressure Drop Through Main Injector Methane Feed')
xlabel('Component')
ylabel('Pressure, bar')
legend('Location','northeast')
% xticks(xs)
xticklabels(Plabels)
grid on
hold off


%% Plot Igniter station conditions
ig_P_idx = [5,7,ig_itr(1:n_ig_pt)];
ig_Ps = logtab{:,ig_P_idx}/1e5;
T_idxs = ig_P_idx+1;
ig_Ts = logtab(:,T_idxs);
xs = 1:length(ig_P_idx);
ys = logtab.t;
[X, Y] = meshgrid(xs,ys);
% labels = logtab.Properties.VariableNames;
Plabels = ["Tank","HVF",ig_sys.PartName{:}];

% figure
% mesh(X,Y,ig_Ps)
% view(45,15)
% title('Pressure Drop Through Igniter Methane Feed')
% xlabel('Component')
% ylabel('Time, s')
% zlabel('Pressure, bar')
% xticks(xs)
% xticklabels(Plabels)

figure
plot(xs,ig_Ps(1,:),'-k','DisplayName','t = 0 s','LineWidth',2)
hold on
plot(xs,interp1(logtab.t,ig_Ps,logtab.t(end)*1/3),'--k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*1/3),'LineWidth',2)
plot(xs,interp1(logtab.t,ig_Ps,logtab.t(end)*2/3),'-.k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)*2/3),'LineWidth',2)
plot(xs,ig_Ps(end,:),':k','DisplayName',sprintf('t = %0.2f s',logtab.t(end)),'LineWidth',2)
title('Pressure Drop Through Igniter Methane Feed')
xlabel('Component')
ylabel('Pressure, bar')
legend('Location','northeast')
% xticks(xs)
xticklabels(Plabels)
grid on
hold off


%% Print key outputs
fprintf('Main Injector mdot: %0.6f kg/s\n',logtab.inj_mdot(1))
fprintf('Main Injector Target mdot: %0.6f kg/s\n',inj_mdot_targ)
fprintf('Main Combustion Chamber Pressure at mdot: %0.2f bar\n',mdot2Pc(logtab.inj_mdot(1))/1e5)
fprintf('Main Injector Outlet Pressure: %0.2f bar\n',logtab.inj_P2(1)/1e5)
fprintf('Main Combustion Chamber Target Pressure: %0.2f bar\n\n',inj_Pc_targ/1e5)
fprintf('Igniter mdot: %0.6f kg/s\n',logtab.ig_mdot(1))
fprintf('Igniter Target mdot: %0.6f kg/s\n',ig_mdot_targ)
fprintf('Igniter Outlet Pressure: %0.2f bar\n',logtab.ig_P2(1)/1e5)
fprintf('Igniter Chamber Target Pressure: %0.2f bar\n',ig_Pc_targ/1e5)
disp(tank_sys)
disp(inj_sys)
disp(ig_sys)