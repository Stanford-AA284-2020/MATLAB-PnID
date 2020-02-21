function Cp = get_Cp(T,i)
% The function returns the Cp fit coefficients for chosen species at T

P = 1e5;

gas = IdealGasMix('gri30.cti');
N     = nSpecies(gas);
iCO  = speciesIndex(gas,'CO');
iH2O = speciesIndex(gas,'H2O');
iCO2  = speciesIndex(gas,'CO2');
iH2   = speciesIndex(gas,'H2');
iCH4   = speciesIndex(gas,'CH4');
xgas = zeros(N,1)';

switch i
    case 1
        
        xgas(iCO) = 1;
        set(gas, 'T', T, 'P', P, 'X', xgas);
        Cp = enthalpy_mole(gas)/1000;
        
    case 2
        
        xgas(iCO2) = 1;
        set(gas, 'T', T, 'P', P, 'X', xgas);
        Cp = enthalpy_mole(gas)/1000;
        
    case 3
        
        xgas(iCH4) = 1;
        set(gas, 'T', T, 'P', P, 'X', xgas);
        Cp = enthalpy_mole(gas)/1000;
        
     case 4
        
        xgas(iH2) = 1;
        set(gas, 'T', T, 'P', P, 'X', xgas);
        Cp = enthalpy_mole(gas)/1000;

     case 5
        
        xgas(iH2O) = 1;
        set(gas, 'T', T, 'P', P, 'X', xgas);
        Cp = enthalpy_mole(gas)/1000;
        
end

end