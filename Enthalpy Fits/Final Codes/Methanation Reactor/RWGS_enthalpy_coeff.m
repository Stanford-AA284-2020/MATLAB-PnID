function x = RWGS_enthalpy_coeff()
%% The function returns coefficients for the best fit to rWGS heat of reaction

gas_h2 = IdealGasMix('gri30.cti');
gas_co = IdealGasMix('gri30.cti');
gas_co2 = IdealGasMix('gri30.cti');
gas_h2o = IdealGasMix('gri30.cti');
gas_ch4 = IdealGasMix('gri30.cti');
P = 1e5;

gas = IdealGasMix('gri30.cti');
N     = nSpecies(gas);
iCO  = speciesIndex(gas,'CO');
iH2O = speciesIndex(gas,'H2O');
iCO2  = speciesIndex(gas,'CO2');
iH2   = speciesIndex(gas,'H2');
iCH4 = speciesIndex(gas,'CH4');
xgas = zeros(1,N);
Range =( (800-100)/0.1 )+ 1;
delta_H = [];
Temp = [];
integrand = [];

%% Constructing the reactions
% CO Methanation


for T = 0:0.1:600
    Ti = T + 273.15;
    i = ((T-0)/0.1) + 1;
    
    xgas = zeros(1,N);
    xgas(iH2) = 1;
    set(gas_h2, 'T', Ti, 'P', P, 'X', xgas);
    h_h2 = enthalpy_mole(gas_h2);
    
    xgas = zeros(1,N);
    xgas(iCO) = 1;
    set(gas_co, 'T', Ti, 'P', P, 'X', xgas);
    h_co = enthalpy_mole(gas_co);
    
    xgas = zeros(1,N);
    xgas(iCO2) = 1;
    set(gas_co2, 'T', Ti, 'P', P, 'X', xgas);
    h_co2 = enthalpy_mole(gas_co2);
    
    xgas = zeros(1,N);
    xgas(iH2O) = 1;
    set(gas_h2o, 'T', Ti, 'P', P, 'X', xgas);
    h_h2o = enthalpy_mole(gas_h2o);
    
    dh = (h_co + h_h2o) - (h_co2 + h_h2);
    delta_H = [delta_H ; dh/1000];
    Temp = [Temp; Ti];
end



F= @(x,xdata)x(1)./(xdata) + x(2).*(xdata.^2) + x(3).*(xdata.^1) + x(4);
x0 = [1,1,1,1];

% F= @(x,xdata)x(1).*(xdata.^(5)) + x(2).*(xdata.^4) + x(3).*(xdata.^3) + x(4).*(xdata.^2)+ x(5).*(xdata.^1) + x(6);
% x0 = [1,1,1,1];
% F= @(x,xdata)x(1).*(xdata.^2) + x(2).*(xdata.^1) + x(3).*(xdata.^0);
% x0 = [1,1,1];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,Temp,delta_H);

figure(4);
T = linspace(273.15, 873.15, 10000);
dh_f = x(1)./(T) + x(2).*(T.^2) + x(3).*(T) + x(4);
plot(T,dh_f, '.k', 'MarkerSize',3);
hold on;
plot(Temp, delta_H, 'b-');
xlabel('Temperature (K)'); ylabel('Heat of reaction (J/mol)');
title('rWGS Enthalpy of Reaction fit');
plotfixer;
end