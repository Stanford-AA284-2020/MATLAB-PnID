% Returns the equilibrium constant for the chosen reaction at T
function Kp = get_eqc(T, i)

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
xgas(iCO2) = 1;
set(gas_co2, 'T', T, 'P', P, 'X', xgas);
g_co2 = gibbs_mole(gas_co2);

xgas = zeros(1,N);
xgas(iH2) = 1;
set(gas_h2, 'T', T, 'P', P, 'X', xgas);
g_h2 = gibbs_mole(gas_h2);

xgas = zeros(1,N);
xgas(iCH4) = 1;
set(gas_ch4, 'T', T, 'P', P, 'X', xgas);
g_ch4 = gibbs_mole(gas_ch4);

xgas = zeros(1,N);
xgas(iCO) = 1;
set(gas_co, 'T', T, 'P', P, 'X', xgas);
g_co = gibbs_mole(gas_co);

xgas = zeros(1,N);
xgas(iH2O) = 1;
set(gas_h2o, 'T', T, 'P', P, 'X', xgas);
g_h2o = gibbs_mole(gas_h2o);

switch i
    
    case 1 %CO Methanation

        delta_G = (g_ch4 + g_h2o) - (g_co + 3*g_h2);
        Kp = exp(-delta_G/(1000*8.314*T));
        
    case 2 %RWGS
        
        delta_G = (g_co + g_h2o) - (g_co2 + g_h2);
        Kp = exp(-delta_G/(1000*8.314*T));
        
    case 3 %CO2 Methanation
        
        delta_G = (g_ch4 + 2*g_h2o) - (g_co2 + 4*g_h2);
        Kp = exp( -delta_G/(1000*8.314*T));
        
end

end