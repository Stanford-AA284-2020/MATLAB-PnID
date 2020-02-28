% Code for computing state evolution of the cryogenic run tank
close all; clear all; clc;
% Provide access to support files via the Matlab path.
addpath 'Fundamental Relation Files' 
addpath 'Fundamental Relation Data'
addpath 'Setup Files' 
addpath 'Property Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n**************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

global Tcrit_i rcrit_i Ttrip_i Tupper_i
global toler
% persistent inumcrit Tnumcrit rnumcrit Pnumcrit

%% Defining the initial conditions %%

% Get the critical and triple data for Oxygen & Nitrogen.
Tc_O2 = Tcrit_i(O2);
rc_O2 = rcrit_i(O2);
Pc_O2 = Pcrit_i(O2);
sc_O2 = s_irT(O2,rc_O2,Tc_O2);
Tt_O2 = Ttrip_i(O2);
Pt_O2 = Ptrip_i(O2);
sft_O2 = s_irT(O2,rftrip_i(O2),Tt_O2);
sgt_O2 = s_irT(O2,rgtrip_i(O2),Tt_O2);

Tc_N2 = Tcrit_i(N2);
rc_N2 = rcrit_i(N2);
Pc_N2 = Pcrit_i(N2);
sc_N2 = s_irT(N2,rc_N2,Tc_N2);
Tt_N2 = Ttrip_i(N2);
Pt_N2 = Ptrip_i(N2);
sft_N2 = s_irT(N2,rftrip_i(N2),Tt_N2);
sgt_N2 = s_irT(N2,rgtrip_i(N2),Tt_N2);

% Fix initial conditions and time step
mdot_0 = 0.21; %kg/s
T0_O2_0 = 140; % K
T0_N2_0 = 300; % K (Assumption: Nitrogen is in equilibrium with the surroundings, and doesn't mix/exchange energy with O2)
P0 = 1000; % psi
P0 = P0*6894.76; % psi to Pa conversion
dt = 0.05; % seconds
V0 = 0.38846; % cubic feet (Internal volume)
V0 = V0*0.0283168; % Convert to cubic meters 

% Compute the initial mass, volume and moles of nitrogen and oxygen
r_N2_0 = rv_iTP(N2, T0_N2_0, P0); % kg/m^2: Nitrogen is expected to be supercritical
r_O2_0 = rl_iTP(O2, T0_O2_0, P0); %  kg/m^3: Oxygen is a subcooled liquid
m_O2_0 = 5; % Fix initial mass of Oxygen (Empirical)
V_O2_0 = m_O2_0/r_O2_0;
V_N2_0 = V0 - V_O2_0;
m_N2_0 = r_N2_0*V_N2_0;  % Initial Mass of Nitrogen
N_N2_0 = m_N2_0/0.028; N_O2_0 = m_O2_0/0.032; % Initial moles of Nirogen and Oxygen

%% Discharge through the tank %%

% Initialize vectors for storing states
time = 0:dt:10; P = zeros(1, length(time)); 
T_O2 =  zeros(1, length(time)); T_N2 =  zeros(1, length(time)); N_O2_L = zeros(1, length(time));
N_O2_V = zeros(1, length(time)); mdot = zeros(1, length(time));downstream_P = zeros(length(time),4); downstream_T = zeros(length(time),4);

T_O2(1,1) = T0_O2_0; T_N2(1,1) = T0_N2_0; P(1,1) = P0; N_O2_L(1,1) = N_O2_0; N_O2_V(1,1) = 0;
m_O2_L = 0.032*N_O2_L; m_O2_V = 0.032*N_O2_V;  [P2_temp, P3_temp, P4_temp, P5_temp, T2_temp, T3_temp, T4_temp, T5_temp] = Downstream_P(O2, P0, mdot_0);
downstream_P(1,:) = [P2_temp, P3_temp, P4_temp, P5_temp]; downstream_T(1,:) = [T2_temp, T3_temp, T4_temp, T5_temp];
mdot(1,1) = mass_flow_rate(O2, P(1,1), T_O2(1,1), downstream_P(1, 1));
i = 1;
for t = dt:dt:10
    t
    i = i + 1;
%     T_O2(1, i-1)
    P_sat = Saturation_iT(O2, T_O2(1, i-1));
    dm = mdot(1, i-1)*dt;
    m_O2_L(1, i) = m_O2_L(1, i-1) - dm;
    r_O2_L = rl_iTP(O2, T_O2(1, i-1), P(1, i-1));
    dq = dm*h_irT(O2, r_O2_L, T_O2(1, i-1)); % Enthalpy carried by the discharged liquid
        
    
    if abs((P(1, i-1) - P_sat)/P_sat) < toler
        % If oxygen is in the liquid-vapor coexistence regime
        disp('Inside the vapordome');
%         v_O2_L = 1/rl_iTP(O2, T_O2(1, i-1), P(1, i-1));
%         
%         Phigh = P(1, i-1); 
%         Plow = 0.8*Phigh;
%         Tnext = Saturation_iP(O2, Phigh);
%         Tlow = Saturation_iP(O2, Plow);
%         Vhigh = (0.028*N_N2_0*rv_iTP(N2, Tnext, Phigh)) + (0.032*n_O2_v*rv_iTP(O2, Tnext, Phigh)) + (0.032*n_O2_l*rl_iTP(O2, Tnext, Phigh));
%         Vlow = (0.028*N_N2_0*rv_iTP(N2, Tlow, Plow)) + (0.032*n_O2_v*rv_iTP(O2, Tlow, Plow)) + (0.032*n_O2_l*rl_iTP(O2, Tlow, Plow));
%         Reshigh = Vhigh - V0;
%         Reslow  = Vlast - V0;
% 
%         if(Reslow*Reshigh > 0)
%             disp('Failed to bracket');
%         end
% 
%         % Start a false position loop to run this to tolerance.
% 
%         NFPs = 200;
%         for j=1:1:NFPs
% 
%             Pnext = Plow + (Phigh - Plow)*(-Reslow)/(Reshigh - Reslow);
%             Tnext = Saturation_iP(O2, Pnext);
% 
%             % Solve Poynting's equation to get mole fraction of oxygen in the
%             % vapor phase
%             x_O2_V = (P_sat/Pnext)*exp((v_O2_L/(8.314*T_O2(1,i-1)/0.032))*(Pnext - P_sat));
%             x_N2_V = 1 - x_O2_V; 
%             N_total = N_N2_0/x_N2_V;
%             n_O2_V = x_O2_V*N_total;
% 
%             % Mass balance to get moles of liquid oxygen remaining
%             n_O2_L = N_O2_L(1, i-1) + N_O2_V(1, i-1) - (dm/0.032) - n_O2_V;
% 
%             % Calculate total volume
%             Vnext = (0.028*N_N2_0*rv_iTP(N2, Tnext, Pnext)) + (0.032*n_O2_v*rv_iTP(O2, Tnext, Pnext)) + (0.032*n_O2_l*rl_iTP(O2, Tnext, Pnext)); 
%             Res   = Vnext - V0;
% 
%             % Check for convergence. Total Volume should be V0
%             if(abs(Res/V0) < toler) 
%                 P_solution = Pnext;          % Found the solution to tolerance.
%                 T_solution = Tnext;
%                 break;
%             end
%             if(Res > 0) % Assumed Pressure is lower
%                 Phigh   = Pnext;    % Midpoint data is on the high side of P.
%                 Reshigh = Res;
%             else
%                 Plow    = Pnext;    % Midpoint data is on the low side of P.
%                 Reslow  = Res;
%             end
%         end
%         P(1, i) = P_solution; T_O2(1, i) = T_solution;
%         N_O2_L(1,i) = n_O2_L; N_O2_V(1,i) = n_O2_V;
%         m_O2_L(1,i) = 0.032*n_O2_L; m_O2_V(1,i) = 0.032*n_O2_V; mdot(1,i) = mdot_0;
        
    else
        
        % Oxygen is still a subcooled liquid
        disp('Subcooled Liquid');
        dV = dm/r_O2_L;
        V_N2_now = (m_N2_0/rv_iTP(N2, T_N2(1, i-1), P(1, i-1))) + dV; % Imposing the volume constraint
        Thigh = T_N2(1, i-1); Tlow = T_N2(1,i-1)*0.8;
        s_N2_high = s_irT(N2, m_N2_0/V_N2_now, Thigh); 
        s_N2_low = s_irT(N2, m_N2_0/V_N2_now, Tlow); 
        s_N2_previous = s_irT(N2, m_N2_0/(V_N2_now-dV), T_N2(1,i-1)); % Entropy at the previous time step
        Reshigh_S = s_N2_high - s_N2_previous;
        Reslow_S = s_N2_low - s_N2_previous;
        
%         if(Reslow_S*Reshigh_S > 0)
%             disp('Failed to bracket');
%         end

        % Start a false position loop to run this to tolerance.

        NFPs = 200;
        for k=1:1:NFPs

            Tnext = Tlow + (Thigh - Tlow)*(-Reslow_S)/(Reshigh_S - Reslow_S);
            snext = s_irT(N2, m_N2_0/V_N2_now, Tnext);
            Res   = snext - s_N2_previous;

            % Check for convergence. Total Volume should be V0
            if(abs(Res/s_N2_previous) < toler) 
                T_solution = Tnext; % Found the solution to tolerance.
                P_solution = P_irT(N2, m_N2_0/V_N2_now, Tnext);
                iteration = k;
                disp('Found Pressure & T_N2');
                break;
            end
            if(Res < 0) % Assumed Temperature is lower. Go higher
                Thigh   = Tnext;    
                Reshigh_S = Res;
            else
                Tlow    = Tnext;    
                Reslow_S  = Res;
            end
        end
        
        P(1, i) = P_solution; T_N2(1, i) = T_solution;
        
        % Update T of LOx now by finding T corresponding to denisty at the
        % updated pressure
        Tlow = T_O2(1, i-1)*0.6; Thigh = T_O2(1, i-1)*1.2;
        r_low = rl_iTP(O2, Tlow, P_solution); r_high = rl_iTP(O2, Thigh, P_solution);
        Reshigh_r = r_high - r_O2_L;
        Reslow_r = r_low - r_O2_L;
        NFPs = 400;
        for l =1:1:NFPs

            Tnext = Tlow + (Thigh - Tlow)*(-Reslow_r)/(Reshigh_r - Reslow_r);
            rnext = rl_iTP(O2, Tnext, P_solution);
            Res   = rnext - r_O2_L;

            % Check for convergence. Total Volume should be V0
            if(abs(Res/r_O2_L) < toler) 
                T_solution_O2 = Tnext; % Found the solution to tolerance.
                iteration = l;
                disp('Found T_LOx');
                break;
            end
            if(Res < 0) % Assumed Temperature is higher. Go lower
                Tlow   = Tnext;    % Midpoint data is on the high side of P.
                Reslow_r = Res;
            else
                Thigh    = Tnext;    % Midpoint data is on the low side of P.
                Reshigh_r  = Res;
            end
        end
       
        if l == NFPs
            T_solution_O2 = T_O2(1, i-1);
        end
    
        m_O2_L(1,i) = m_O2_L(1, i-1) - dm; m_O2_V(1,i) = 0; 
        N_O2_L(1,i) = m_O2_L(1,i)/0.032; N_O2_V(1,i) = 0; T_O2(1, i) = T_solution_O2;
        mass_flow_rate(O2, P(1,i), T_O2(1,i), downstream_P(i-1, 1))
        mdot(1,i) = mass_flow_rate(O2, P(1,i), T_O2(1,i), downstream_P(i-1, 1));
        [P2_temp, P3_temp, P4_temp, P5_temp, T2_temp, T3_temp, T4_temp, T5_temp] = Downstream_P(O2, P(1,i), mdot(1,i));
        downstream_P(i,:) = [P2_temp, P3_temp, P4_temp, P5_temp];
        downstream_T(i,:) = [T2_temp, T3_temp, T4_temp, T5_temp];

    end
      
  
end

save('results.mat', 'm_O2_L', 'T_O2', 'mdot', 'downstream_P', 'downstream_T', 'P');

figure(1);
clf;
plot(time, P/1e5, 'k-');
xlabel('Time [s]'); ylabel('Pressure [bar]');
title('Run tank pressure during blowdown');
axis([0 10 20 80]);
plotfixer;

figure(2);
clf;
plot(time, m_O2_L, 'k-');
xlabel('Time [s]'); ylabel('Mass [kg]');
title('Mass of LOX remaining during blowdown');
axis([0 10 0 8]);
plotfixer;

figure(3);
clf;
plot(time, mdot, 'k-');
xlabel('Time [s]'); ylabel('Mass flow rate [kg/s]');
title('Mass flow rate through the LOX line');
axis([0 10 0 0.3]);
plotfixer;

figure(4);
clf;
plot(time, T_O2, 'k-');
xlabel('Time [s]'); ylabel('Temperature [K]');
title('Temperature of LOX in tank');
axis([0 10 120 160]);
plotfixer;

figure(5);
clf;
hold on;
t1_ind = 51; t2_ind = 101; t3_ind = 151; t4_ind = 201;
P1_vect = [P(1, t1_ind), downstream_P(t1_ind, :)]/1E5;
P2_vect = [P(1, t2_ind), downstream_P(t2_ind, :)]/1E5;
P3_vect = [P(1, t3_ind), downstream_P(t3_ind, :)]/1E5;
P4_vect = [P(1, t4_ind), downstream_P(t4_ind, :)]/1E5;
x_stations = [1, 2, 3, 4, 5];
axis([1, 5, 0, 80]);
plot(x_stations, P1_vect, 'k-');
plot(x_stations, P2_vect, 'k--');
plot(x_stations, P3_vect, 'k.-');
plot(x_stations, P4_vect, 'k:');
legend('t = 2.5 s', 't = 5 s', 't = 7.5 s', 't = 10 s');
ax = gca;
ax.XTick = x_stations;
ax.XTickLabels = {'Tank', 'CVMO', 'CKMO', 'BVMO', 'Injector'};
xlabel('Station'); ylabel('Pressure [bar]');
title('Pressure ladder - LOx Feed Line');
plotfixer;
