function output = regulator(medium, P1, T1, Cv, P2reg)
% regulator Outlet P & T
%   output = regulator(medium, P1, T1, Cv, P2reg) Calculates pressure and
%   temperature drop across a flow control device defined by its flow
%   coefficient Cv and regulated outlet pressure P2reg. If regulator is 
%   choked, returns choked mass flow rate instead of pressure and 
%   temperature downstream.
%
%   Formulas for dP based on Cv are from Swagelok Technical
%   Bulletin MS-06-84-E Rev. 4 (2007), which is itself derived from ISA 
%   S75.01, Flow Equations for Sizing Control Valves, Standards and 
%   Recommended Practices for Instrumentation and Control, 10th ed., 
%   Vol. 2, 1989.

    rho1 = PREoS(medium, "rho", P1, T1);% Fluid density at P1, T1
    
    if length(rho1) > 1
        error('Mixed Phases - Check for Cavitation')
    elseif rho1 > 500 % If likely liquid
        error('Liquid flow through regulator')
        
    else
        % Calculate pressure drop across valve for compressible fluid
        N2 = 6950;% Unit coefficient, bar & std. L/min
        P1b = P1/1e5;% bar, Converted from Pa
        Gg = medium.Mw/28.96;% Specific Gravity relative to Air
        SLPM2kgps = PREoS(medium,"rho",101325,273.15)/(1000*60);
        q = mdot/SLPM2kgps;% SLPM Vol flow
        q_choke = 0.471*N2*Cv*P1b*sqrt(1/(Gg*T1));% Flow thru choked valve
        if q > q_choke% If valve flow is choked, output choked mass flow
            output = q_choke*SLPM2kgps;
            return
        end
        
        % If not choked, P2 = P2reg
        % ASSUMING regulator is actively controlled to eliminate droop
        P2 = P2reg;
        % Calculate outlet temperature, assuming isentropic expansion
        T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
        output = [P2, T2];
        
    end
end