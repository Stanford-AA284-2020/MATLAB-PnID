function x = CO_Methanation_enthalpy_coeff2()
%% The function returns coefficients for the best fit to CO Methanation's heat of reaction

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


for T = 100:0.1:800
    Ti = T + 273.15;
    i = ((T-200)/0.1) + 1;
    
    xgas = zeros(1,N);
    xgas(iH2) = 1;
    set(gas_h2, 'T', Ti, 'P', P, 'X', xgas);
    h_h2 = enthalpy_mole(gas_h2);
    
    xgas = zeros(1,N);
    xgas(iCH4) = 1;
    set(gas_ch4, 'T', Ti, 'P', P, 'X', xgas);
    h_ch4 = enthalpy_mole(gas_ch4);
    
    xgas = zeros(1,N);
    xgas(iCO) = 1;
    set(gas_co, 'T', Ti, 'P', P, 'X', xgas);
    h_co = enthalpy_mole(gas_co);
    
    xgas = zeros(1,N);
    xgas(iH2O) = 1;
    set(gas_h2o, 'T', Ti, 'P', P, 'X', xgas);
    h_h2o = enthalpy_mole(gas_h2o);
    
    dh = -(h_ch4 + h_h2o) + (h_co + 3*h_h2);
    delta_H = [delta_H ; dh/1000];
    Temp = [Temp; Ti];
end

plot(Temp, delta_H, 'b-')
hold on;

F= @(x,xdata)x(1).*(xdata.^2) + x(2).*(xdata.^1) + x(3).*(xdata.^0);
x0 = [1,1,1];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,Temp,delta_H);

end