% Main Injector Performance

close all; clear all;

addpath(genpath('/Users/JBR132/Documents/_Stanford/AA284B Propulsion System Design Lab/MATLAB-PnID'))

load('methane_feed_system_output.mat')
load('results.mat')

FluidDatabase

%% Fluid Conditions

T_Ox = downstream_T(:,3);% Saturation pressure at 10 bar
P_Ox = downstream_P(:,3);
T_CH4 = logtab.BVMF2_T2(:);% Gas temp upstream of injector
P_CH4 = logtab.BVMF2_P2(:);

mdot_Ox = mdot;% kg/s
mdot_CH4 = logtab.inj_mdot(:);% kg/s

Pc = 10e5;% Pa


%% Injector Geometry

% A_CH4 = inj_sys.A(end);% Get directly from system
A_Ox = orifice_size(Oxygen,mdot_Ox(1),mean(P_Ox(1:20)) - Pc,0.6,mean(P_Ox(1:20)),mean(T_Ox(1:20)))