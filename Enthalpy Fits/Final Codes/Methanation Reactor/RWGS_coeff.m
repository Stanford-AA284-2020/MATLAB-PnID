function x = RWGS_coeff()
%% The function returns coefficients for the best fit to RWGS equilbirum constant

gas_h2 = IdealGasMix('gri30.cti');
gas_co = IdealGasMix('gri30.cti');
gas_co2 = IdealGasMix('gri30.cti');
gas_h2o = IdealGasMix('gri30.cti');
P = 1e5;

gas = IdealGasMix('gri30.cti');
N     = nSpecies(gas);
iCO  = speciesIndex(gas,'CO');
iH2O = speciesIndex(gas,'H2O');
iCO2  = speciesIndex(gas,'CO2');
iH2   = speciesIndex(gas,'H2');
xgas = zeros(1,N);
Range =( (800-100)/0.1 )+ 1;
Kp = [];
Temp = [];
integrand = [];

%% Constructing the reactions
% Reverse Water Gas Shift Reaction


for T = 0:0.1:500
    Ti = T + 273.15;
    i = ((T-0)/0.1) + 1;
    
    xgas = zeros(1,N);
    xgas(iH2) = 1;
    set(gas_h2, 'T', Ti, 'P', P, 'X', xgas);
    g_h2 = gibbs_mole(gas_h2);
    
    xgas = zeros(1,N);
    xgas(iCO2) = 1;
    set(gas_co2, 'T', Ti, 'P', P, 'X', xgas);
    g_co2 = gibbs_mole(gas_co2);
    
    xgas = zeros(1,N);
    xgas(iCO) = 1;
    set(gas_co, 'T', Ti, 'P', P, 'X', xgas);
    g_co = gibbs_mole(gas_co);
    
    xgas = zeros(1,N);
    xgas(iH2O) = 1;
    set(gas_h2o, 'T', Ti, 'P', P, 'X', xgas);
    g_h2o = gibbs_mole(gas_h2o);
    
    delta_G = -(g_co2 + g_h2) + (g_co + g_h2o);
    kp = -delta_G/(8.314*Ti);
    Kp = [Kp ; exp(kp/1000)];
    Temp = [Temp; Ti];
end

log_Kp = log(Kp);


F = @(x,xdata)x(1)./(xdata.^2) + x(2)./(xdata) + x(3);
x0 = [1,1,1];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,Temp,log_Kp);
figure(1);
T = linspace(273.15, 773.15, 10000);
Kp_f = exp(x(1)./(T.^2) + x(2)./(T) + x(3));

figure(3);
plot(T,Kp_f, '.k', 'MarkerSize',3);
hold on;
plot(Temp, Kp, 'b-');
xlabel('Temperature (K)'); ylabel('Kp');
title('rWGS Kp fit');

end