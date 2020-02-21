%% This function constructs the species conservation equations

function dydt = methanation_mk9(t,y)

%T = 70 + 273.14; %Reactor temperature
global velocity 
global rho_cat 
global P
global Q
global Ru


%Check for negative temperature or negative hydrogen concentration

% if y(1) < 0 || y(5) < 0
%     disp('Negative y(1)erature or Concentration: Returning');
%     return;
% end

%Concentration to Mole Fraction
z_CO = y(2)/(P/(Ru*y(1)));
z_CO2 = y(3)/(P/(Ru*y(1)));
z_CH4 = y(4)/(P/(Ru*(y(1))));
z_H2 = y(5)/(P/(Ru*(y(1))));
z_H2O = y(6)/(P/(Ru*y(1)));
mole_sum = z_CO + z_CO2 + z_CH4 + z_H2 + z_H2O;
z_CO = z_CO/mole_sum;
z_CO2 = z_CO2/mole_sum;
z_CH4 = z_CH4/mole_sum;
z_H2 = z_H2/mole_sum;
z_H2O = z_H2O/mole_sum;

load('best_fits7.mat');
x1 = x(1, :); x2 = x(2, :); x3 = x(3, :); 
x4 = x(4, :); x5 = x(5, :);
%get_Cp should return the Cp values in J/mol-K
Cp_H2 = x1(1).*(y(1).^2) + x1(2).*(y(1).^1) + x1(3);
Cp_CO2 = x2(1).*(y(1).^2) + x2(2).*(y(1).^1) + x2(3);
Cp_CO = x3(1).*(y(1).^2) + x3(2).*(y(1).^1) + x3(3);
Cp_H2O = x4(1).*(y(1).^2) + x4(2).*(y(1).^1) + x4(3);
Cp_CH4 = x5(1).*(y(1).^2) + x5(2).*(y(1).^1) + x5(3);


Cp_avg = ((z_H2*Cp_H2) + (z_CO2*Cp_CO2) + (z_CO*Cp_CO) + (z_CH4*Cp_CH4) + (z_H2O*Cp_H2O)); % J/mol-K

%Total molecular weight
M = (28*z_CO + 44*z_CO2 + 16*z_CH4 + 2*z_H2 + 18*z_H2O)/1000;
Cp_avg = Cp_avg/M; % J/kg-K

%% Description of rate constants
K_CO = (3.14*(10^(-15)))*exp(-120200/(Ru*y(1))); %Adsorption constant for CO, 1/Pa
K_H2 = (6.12*(10^(-4)))*exp(-82900/(Ru*y(1))); %Adsorption constant for H2, 1/bar
K_CH4 = (6.65*(10^(-9)))*exp(-38280/(Ru*y(1))); %Adsorption constant for CH4, 1/bar
K_H2O = (1.77*(10^(5)))*exp(88680/(Ru*y(1))); %Dissociative adsorption constant of H2, 1/bar??

k1 = (1.83*(10^(15)))*exp(-218900/(Ru*y(1))); %Rate coefficient of Reaction 1, kmol.(bar^0.5)/(kg-cat. hr)
k2 = (31.1)*exp(-59400/(Ru*y(1))); %Rate coefficient of Reaction 2, kmol/(kg-cat. hr)59400
k3 = (4.05*(10^(15)))*exp(-209900/(Ru*y(1))); %Rate coefficient of Reaction 3, kmol.(bar^0.5)/(kg-cat. hr)

K1 = 0.0005*exp((x_K1(1)/((y(1))^2)) + (x_K1(2)/((y(1))^1)) + x_K1(3))*((P)^2); %Equilibrium constant of Reaction 1 (Carbon Monoxide Methanation), bar^2, from Zhang
K2 = 500*exp((x_K2(1)/((y(1))^2)) + (x_K2(2)/((y(1))^1)) + x_K2(3)); %Equilibrium constant of Reaction 2 (Water Gas Shift)
K3 = 0.8*exp((x_K3(1)/((y(1))^2)) + (x_K3(2)/((y(1))^1)) + x_K3(3))*((P)^2); %Equilibrium constant of Reaction 1 (Carbon Dioxide Methanation), bar^2



DEN = 1 + (K_CO*z_CO*P) + (K_H2*z_H2*P) + (K_CH4*z_CH4*P)+ (K_H2O*z_H2O/z_H2);

r1 = -(k1/((z_H2*P)^2.5))*(((z_CH4*P)*(z_H2O*P)) - (((z_H2*P)^3)*(z_CO*P)*K1))/(DEN^2);
r2 = -(k2/((z_H2*P)^1))*(((z_CO*P)*(z_H2O*P)) - (((z_H2*P)^1)*(z_CO2*P)*K2))/(DEN^2);
r3 = -(k3/((z_H2*P)^3.5))*(((z_CH4*P)*((z_H2O*P)^2)) - (((z_H2*P)^4)*(z_CO2*P)*K3))/(DEN^2);

% delta_H1 = get_hr(y(1), 1);
% delta_H2 = get_hr(y(1), 2);
% delta_H3 = get_hr(y(1), 3); % J/mol
delta_H1 = (x_H1(1)*((y(1))^(-1))) +  (x_H1(2)*((y(1))^2)) + (x_H1(3)*((y(1))^1)) + (x_H1(4))
delta_H2 = (x_H2(1)*((y(1))^(-1))) +  (x_H2(2)*((y(1))^2)) + (x_H2(3)*((y(1))^1)) + (x_H2(4))
delta_H3 = (x_H3(1)*((y(1))^(-1))) +  (x_H3(2)*((y(1))^2)) + (x_H3(3)*((y(1))^1)) + (x_H3(4))

% m_flux = rho*velocity;
c_tot = (z_CO + z_CO2 + z_H2 + z_H2O + z_CH4)*P/(Ru*y(1));
rho = P/(c_tot*Ru*y(1)/M);

%%Setting up the differential equations
dydt = zeros(6,1);
dydt(1) = rho_cat*((r1*(-delta_H1)) + (r2*(-delta_H2)) + (r3*(-delta_H3)))/(Cp_avg*velocity*c_tot) - (Q/(rho*Cp_avg*velocity));% Convective Heat Transfer neglected
dydt(2) = rho_cat*(- r1 + r2 )/velocity;
dydt(3) = rho_cat*(-r3 - r2)/velocity;
dydt(4) = rho_cat*(r1 + r3)/velocity;
dydt(5) = -rho_cat*(3*r1 + r2 + 4*r3)/velocity;
dydt(6) = -rho_cat*(-r1 - r2 - 2*r3)/velocity

end
