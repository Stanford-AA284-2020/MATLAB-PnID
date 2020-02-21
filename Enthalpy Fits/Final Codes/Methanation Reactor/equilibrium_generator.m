gas = IdealGasMix('gri30.cti');
N     = nSpecies(gas);
iCO  = speciesIndex(gas,'CO');
iH2O = speciesIndex(gas,'H2O');
iCO2  = speciesIndex(gas,'CO2');
iH2   = speciesIndex(gas,'H2');
iCH4 = speciesIndex(gas,'CH4');
xgas = zeros(1,N);
CO = []; CO2 =[]; CH4 = []; H2O = []; H2 = []; Tf = [];
%Equilibrium for Methanation
for T1 = 273:1:873; 
P1 = 2e7;
xgas(iCO) = 0.2;
xgas(iCO2) = 0; 
xgas(iH2) = 0.8;
xgas(iH2O) = 0;
gas1 = gas;
set(gas1, 'T', T1, 'P', P1, 'X', xgas);
h1 = enthalpy_mole(gas1)/1000; 

gas2 = equilibrate(gas1, 'HP');
xgas2 = moleFractions(gas2);
CO_eq = xgas2(iCO);
CO2_eq = xgas2(iCO2);
CH4_eq = xgas2(iCH4);
H2_eq = xgas2(iH2);
H2O_eq = xgas2(iH2O);

CO = [CO; CO_eq];
CO2 = [CO2; CO2_eq];
CH4 = [CH4; CH4_eq];
H2 = [H2; H2_eq];
H2O = [H2O; H2O_eq];
Tf = [Tf; T1]
end

figure(1);
hold on;
plot(Tf, CO, '-g');
plot(Tf, CO2, '-y');
plot(Tf, CH4, '-r');
plot(Tf, H2, '-k');
plot(Tf, H2O, '-b');
%text(0.15, 0.32, {' \it Dotted Lines : Equilibrium concentrations '});
xlabel('Temperature (K)');
ylabel('Mole Fractions');
legend('CO', 'CO2', 'CH4', 'H2', 'H2O');
plotfixer;
