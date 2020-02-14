function output = valve(medium, mdot, P1, T1, Cv)
% valve Outlet P & T
%   [P2, T2] = valve(medium, mdot, P1, T1, Cv) Calculates pressure and
%   temperature drop across a flow control device defined by its flow
%   coefficient Cv. If valve is choked, returns choked flow rate instead of
%   pressure and temperature downstream.
%
%   Formulas for dP based on Cv are from Swagelok Technical
%   Bulletin MS-06-84-E Rev. 4 (2007), which is itself derived from ISA 
%   S75.01, Flow Equations for Sizing Control Valves, Standards and 
%   Recommended Practices for Instrumentation and Control, 10th ed., 
%   Vol. 2, 1989.

    rho1 = PREoS(medium, "rho", P1, T1);% Fluid density at P1, T1
    if length(rho1) > 1
        error('Mixed Phases - Check for Cavitation')
    end
    
    if rho1 > 500% If likely liquid, use liquid Cv formula from Swagelok
        % Calculate pressure drop across valve for incompressible fluid
        N1 = 14.42;% Unit coefficient, bar & L/min
        Gf = rho1/1000;% Specific Gravity relative to Water at 4C
        q = mdot/rho1*1000*60;% Volumetric flow rate, L/min
        % q = N1*Cv*sqrt(dP/Gf);% Formulation in Swagelok Tech Bulletin
        dP = Gf*(q/(N1*Cv))^2;% bar, Rearranged to output dP
        dP = dP*1e5;% Pa, Converted from bar
        P2 = P1-dP;
        T2 = PREoS(medium, "T", P2, rho1);% Incompressible flow assumption?
        output = [P2, T2];
        
    else
        % Calculate pressure drop across valve for compressible fluid
        N2 = 6950;% Unit coefficient, bar & std. L/min
        P1b = P1/1e5;% bar, Converted from Pa
        Gg = medium.Mw/28.96;% Specific Gravity relative to Air
        SLPM2mdot = PREoS(medium,"rho",101325,273.15)/(1000*60);
        q = mdot/SLPM2mdot;% SLPM Vol flow
        q_choke = 0.471*N2*Cv*P1b*sqrt(1/(Gg*T1));% Flow thru choked valve
        if q > q_choke% If valve flow is choked, output choked mass flow
            output = q_choke*SLPM2mdot;
            return
        end
        
        % Find dP numerically within non-choked bounds
        % Objective Function: Formulation in Swagelok Tech Bulletin - q
        qfn = @(dP) abs(N2*Cv*P1b*(1-2*dP/(3*P1b))*sqrt(dP/(P1b*Gg*T1))-q);
        dP = fminbnd(qfn,0,P1b*(1-1/1.8));% Golden section search for min
        P2 = P1-dP*1e5;% Pa, Converted from bar
        % Calculate outlet temperature, assuming isentropic expansion
        T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
        output = [P2, T2];
        
    end
end