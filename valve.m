function [P2, T2] = valve(medium, mdot, P1, T1, Cv)
    rho1 = PREoS(medium, "rho", P1, T1);% Fluid density at P1, T1
    if length(rho1) > 1
        error('Mixed Phases Detected, Check for Cavitation')
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
        T2 = T1;% Incompressible flow assumption?
    else
        % Calculate pressure drop across valve for compressible fluid
        N2 = 6950;% Unit coefficient, bar & std. L/min
        P1b = P1/1e5;% bar, Converted from Pa
        Gg = medium.Mw/28.96;% Specific Gravity relative to Air
        q = mdot/PREoS(medium,"rho",101325,273.15)*1000*60;% SLPM Volumetric flow rate
        qfn = @(dP) abs(N2*Cv*P1b*(1-2*dP/(3*P1b))*sqrt(dP/(P1b*Gg*T1)) - q);% Formulation in Swagelok Tech Bulletin
        dP = fminbnd(qfn,0,P1b/1.8);% Find dP numerically within non-choked bounds
        P2 = P1-dP*1e5;% Pa, Converted from bar
        % Calculate outlet temperature, assuming isentropic expansion
        T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
    end
end