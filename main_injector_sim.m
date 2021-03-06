% Main Injector Performance

close all; clear all;

addpath(genpath('/Users/JBR132/Documents/_Stanford/AA284B Propulsion System Design Lab/MATLAB-PnID'))

load('methane_feed_system_output.mat')
load('results.mat')

FluidDatabase


%% Fluid Conditions

T_Ox = downstream_T(:,3);% Saturation pressure at 10 bar
P_Ox = downstream_P(:,3);
T_CH4 = logtab.eleminlets_T2(:);% Gas temp upstream of injector
P_CH4 = logtab.eleminlets_P2(:);

% mdot_Ox = mdot';% kg/s
mdot_CH4 = logtab.inj_mdot(:);% kg/s

Pc = 10e5;% Pa


%% Injector Geometry
n_elem = 5;
minwall = 0.015*0.0254;% in -> m, from Protolabs design guide min feature size
mingap = 0.025*0.0254;% in -> m, from Protolabs feedback

% A_Ox = orifice_size(Oxygen,mdot_Ox(1),mean(P_Ox(1:20)) - Pc,0.6,mean(P_Ox(1:20)),mean(T_Ox(1:20)))
% A_Ox = 1e-5;
A_Ox = pi*(0.004^2)/4;% From Downstream_P function in Cryo Code Mk7
A_Oxelem = A_Ox/n_elem;
D_Oxelem = sqrt(A_Oxelem/pi)*2% m

ID_CH4elem = D_Oxelem + 2*minwall
OD_CH4elem = ID_CH4elem + 2*mingap*1.2% extra margin on min gap to avoid plugging
A_CH4elem = pi*(OD_CH4elem^2 - ID_CH4elem^2)/4;
A_CH4 = A_CH4elem*n_elem
% A_CH4 = inj_sys.A(end);% Get directly from system
% A_CH4elem = A_CH4/n_elem;

rho_Ox = zeros(length(logtab.t),1);
rho_CH4 = zeros(length(logtab.t),1);
mdot_Ox = zeros(length(logtab.t),1);
tb = zeros(length(logtab.t),1);
% V_Ox = zeros(length(logtab.t),1);
% V_CH4 = zeros(length(logtab.t),1);
t_Ox = [0:0.05:10]';% time steps for ox sim
t_CH4 = logtab.t;
for step = 1:length(t_CH4)
    t = t_CH4(step);
    rho_Ox(step) = PREoS(Oxygen,"rho",interp1(t_Ox,P_Ox,t),interp1(t_Ox,T_Ox,t));
    rho_CH4(step) = PREoS(Methane,"rho",interp1(t_CH4,P_CH4,t),interp1(t_CH4,T_CH4,t));
    mdot_Ox(step) = interp1(t_Ox,mdot,t);
%     V_Ox(step) = mdot_Ox(step)/(rho_Ox(step)*A_Ox);
%     V_CH4(step) = mdot_CH4(step)/(rho_CH4(step)*A_CH4);
%     tb(step) = 5*O2DynVisc(interp1(t_Ox,T_Ox,t),rho_Ox(step))/(rho_Ox(step)*(V_CH4(step)-V_Ox(step))^2);
%     [P2_CH4, T2_CH4] = 
end

V_Ox = mdot_Ox./(rho_Ox*A_Ox);
V_CH4 = mdot_CH4./(rho_CH4*A_CH4);

VR = V_CH4./V_Ox;
J = (rho_CH4.*V_CH4.^2)./(rho_Ox.*V_Ox.^2);

% tb = (D_Oxelem/2)./(V_CH4 - V_Ox) .* sqrt(rho_Ox./rho_CH4);
% Lb = tb.*V_CH4;
Lb = 25*D_Oxelem./J.^0.2;

subplot(1,3,1)
plot(t_CH4,VR,'-k','LineWidth',2,'DisplayName','Velocity Ratio')
title('Velocity Ratio')
xlabel('Time, s')
subplot(1,3,2)
plot(t_CH4,J,'-k','LineWidth',2,'DisplayName','Momentum Flux Ratio')
title('Momentum Flux Ratio')
xlabel('Time, s')
subplot(1,3,3)
plot(t_CH4,Lb,'-k','LineWidth',2,'DisplayName','Jet Breakup Distance')
title('Jet Breakup Distance (m)')
xlabel('Time, s')

function n = O2DynVisc(T,rhoi)
    rho = rhoi/1e3;% kg/m3 -> g/cm3
    consts = [8.3447128902e-7;...
        -6.4319584704e-9;...
        1.3385840412e-10;...
        -1.3283970709e-12;...
        7.4604203289e-15;...
        -2.5398419748e-17;...
        5.2077731498e-20;...
        -5.9258985472e-23;...
        2.8773126492e-26];
    n0 = 0;
    for i=1:length(consts)
        n0 = n0 + consts(i)*T^i;
    end
    if rho > 0.932 % g/cm3
        nE = (0.65391907848*rho + 0.000029886313449*(exp(9.25*rho)-1))*10^-3;% g/cm-s
    else
        nE = 4.7293376329e-4*rho + -1.7410413651e-4*rho^2 + 5.9995361171e-4*rho^3;% g/cm-s
    end
    n = (n0 + nE)*100/1000; % g/cm-s -> kg/m-s
end