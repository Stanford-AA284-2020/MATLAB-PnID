function output = minorloss(medium, mdot, P1, T1, K, A, temp_interp)
% minorloss Pressure and temperature downstream of minor loss
%   output = minorloss(medium, mdot, P1, T1, K, A, temp_interp) Calculates 
%   pressure and temperature drop across a minor hydraulic loss 
%   characterized by a minor loss coefficient K and orifice area A. If
%   choked, output is empty since these equations assume incompressible
%   flow.
%
%   temp_interp = "isothermal", "adiabatic", "isentropic" specifies the
%   method used to find the temperature downstream of the orifice given the
%   pressure drop
%
%   Formula for dP from Cd, A from Fox & McDonald's Introduction to Fluid
%   Mechanics (8e) Section 8.10: Restriction Flow Meters for Internal 
%   Flows, pp. 389 (beta^4 in denominator assumed negligibly small)

    rho1 = PREoS(medium, "rho", P1, T1);% Fluid density at P1, T1
    
    if length(rho1) > 1
        error('Mixed Phases - Check for Cavitation')
        
    else
        v = mdot/(rho1*A);% Flow velocity at inlet of minor loss
        P2 = P1 - K/2*rho1*v^2; % Fox & McDonald 8.40a, Feedline 284b slide
        if P2 <= P1/chokeratio(medium.gam)
            output = [];
            return
        end
        % Calculate outlet temp using desired method
        T2 = T2_interp(medium, P1, T1, P2, temp_interp);
        output = [P2, T2];
        
    end
end