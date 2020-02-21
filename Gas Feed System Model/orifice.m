function output = orifice(medium, mdot, P1, T1, Cd, A, temp_interp)
% orifice Outlet P & T
%   output = orifice(medium, mdot, P1, T1, Cd, A, temp_interp) Calculates 
%   pressure and temperature drop across a flow control device defined by 
%   its discharge coefficient Cd and orifice area A. If orifice is choked, 
%   returns choked mass flow rate instead of pressure and temperature.
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
        
    elseif rho1 > 500 % If likely liquid
        error('Liquid flow through regulator')
        
    else
        % Calculate pressure drop across orifice for compressible fluid
        R = 8314.46261815324;% J/kmol-K, Gas Constant
        k = medium.gam;% Specific Heat Ratio
        Mw = medium.Mw;% kg/kmol, Molar Mass
        mdot_choke = k/((k+1)/2)^((k+1)/(2*(k-1)))*(P1*A/sqrt(k*R/Mw*T1));
        if mdot > mdot_choke % If orifice is choked, output choked mass flow
            output = mdot_choke;
            return
        end
        
        % If not choked, calculate P2 from Fox & McDonald Eq. 8.54
        % mdot = Cd*A*sqrt(2*rho1*(P1-P2))
        P2 = P1 - (mdot/(Cd*A))^2/(2*rho1);
        % Calculate outlet temperature, assuming isentropic expansion
%         T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
        % Calculate outlet temp using desired method
        T2 = T2_interp(medium, P1, T1, P2, temp_interp);
        output = [P2, T2];
        
    end
end