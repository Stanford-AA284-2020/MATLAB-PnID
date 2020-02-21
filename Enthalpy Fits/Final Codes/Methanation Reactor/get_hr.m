function delta_h = get_hr(T, i)
% The function returns the heat of reaction for the chosen reaction at T

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

xgas = zeros(1,N);
xgas(iCO2) = 1;
set(gas_co2, 'T', T, 'P', P, 'X', xgas);
h_co2 = enthalpy_mole(gas_co2);
  
xgas = zeros(1,N);
xgas(iH2) = 1;
set(gas_h2, 'T', T, 'P', P, 'X', xgas);
h_h2 = enthalpy_mole(gas_h2);

xgas = zeros(1,N);
xgas(iCH4) = 1;
set(gas_ch4, 'T', T, 'P', P, 'X', xgas);
h_ch4 = enthalpy_mole(gas_ch4);

xgas = zeros(1,N);
xgas(iCO) = 1;
set(gas_co, 'T', T, 'P', P, 'X', xgas);
h_co = enthalpy_mole(gas_co);

xgas = zeros(1,N);
xgas(iH2O) = 1;
set(gas_h2o, 'T', T, 'P', P, 'X', xgas);
h_h2o = enthalpy_mole(gas_h2o);

switch i
    
    case 1 % CO Methanation
        
        dh = (h_ch4 + h_h2o) - (h_co + 3*h_h2);
        delta_h = dh/1000;

    case 2 % RWGS
        
        dh = (h_co + h_h2o) - (h_co2 + h_h2);
        delta_h = dh/1000;
           
    case 3 % CO2 Methanation
        
        dh = (h_ch4 + 2*h_h2o) - (h_co2 + 4*h_h2);
        delta_h = dh/1000;
        
end

end

