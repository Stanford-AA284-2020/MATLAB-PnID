function out = PREoS(medium, find, in1, in2)
% PREoS Peng-Robinson Equation of State for Pure Fluids
%   out = PREoS(medium, find, in1, in2) Calculates the remainder of
%   Pressure, Temperature, and Density depending on the selected inputs.
%   
%   medium: Methane, Oxygen, Nitrogen, Air, Water
%   Input variables according to following table:
%   find | "rho" |  "P"  |  "T"  |
%    in1 |  "P"  |  "T"  |  "P"  |
%    in2 |  "T"  | "rho" | "rho" |
%   P: Pressure, Pa
%   T: Temperature, K
%   rho: Density, kg/m^3
%
%   Sources:
%   D. Peng and D. Robinson, "A New Two-Constant Equation of State", 
%       Ind. Eng. Chem.,Fundam., Vol. 15, No. 1, 1976, pp. 59-64
%   Section 2: Physical and Chemical Data, Perry's Chemical Engineers' 
%       Handbook, Seventh Edition

% Gas Constant
R = 8314.46261815324;% J/kmol-K

% Get Fluid Properties from Medium
Mw = medium.Mw;
Tc = medium.Tc;
Pc = medium.Pc;
w = medium.w;

% Calculate parameters applicable to all outputs
a_Tc = 0.45724*R^2*Tc^2/Pc;% a at critical temp
b = 0.07780*R*Tc/Pc;% b at critical temp = b at desired temp
kap = 0.37464 + 1.54226*w - 0.26992*w^2;% Characteristic constant

% Calculate Output Requested
if find == "P"
    T = in1;% K
    rho = in2;% kg/m^3
    Vm = Mw/rho;% m^3/kmol
    % Use typical form (Eq 4) of PREoS to find P explicitly
    P = preos(T, Vm, b);
    out = P;
    
elseif find == "rho"
    P = in1;% Pa
    T = in2;% K
    A = a(T)*P/(R^2*T^2);
    B = b*P/(R*T);
    % Use cubic form of PREoS to find compressibility factor
    Z = roots([1; B-1; A-3*B^2-2*B; B^3+B^2-A*B]);
    Z = Z(imag(Z)==0); % Discard imaginary roots
    if length(Z) == 1% If gas only
        Vm = Z(1)*R*T/P;% m^3/kmol
        rho = Mw/Vm;% kg/m^3
        out = rho;
        
    elseif length(Z) > 1% If vapor & liquid
        % Vapor
        Vm_vap = max(Z)*R*T/P;% m^3/kmol
        rho_vap = Mw/Vm_vap;% kg/m^3
        % Liquid
        Vm_liq = min(Z)*R*T/P;% m^3/kmol
        rho_liq = Mw/Vm_liq;% kg/m^3
        out = [rho_vap; rho_liq];
        
    end
    
elseif find == "T"
    P = in1;% Pa
    rho = in2;% kg/m^3
    Vm = Mw/rho;% m^3/kmol
    % Find roots of typical form (Eq 4) of PREoS for T
    zero = @(T) preos(T, Vm, b) - P;
    T = fzero(zero,Tc);
    out = T;
    
end

% a coefficient (used in all 3 
function a_T = a(T)
    alpha = (1 + kap*(1-sqrt(T/Tc)))^2;% Dimensionless fxn of reduced temp
    a_T = a_Tc*alpha;% a at desired temp
end

% Conventional form of Peng-Robinson EoS (Peng Eq. 4)
function P = preos(T, Vm, b)
    P = (R*T)/(Vm - b) - a(T)/(Vm*(Vm+b) + b*(Vm - b));
end
end