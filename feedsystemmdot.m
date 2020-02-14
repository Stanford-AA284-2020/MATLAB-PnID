function output = feedsystemmdot(systable, medium, Ptank, Ttank, A, mdot, varargin)
% Propellant Feed System as function of mass flow
%   Takes system definition table, working fluid, tank conditions, control
%   orifice area, and mass flow rate.
%   Calculates fluid conditions downstream of each component up to orifice,
%   then calculates mass flow that would correspond to fluid conditions
%   just upstream of orifice.
%   Returns error between input mass flow rate and resultant orifice flow
%   rate, the orifice flow rate itself, and the updated system table for
%   the given input mass flow rate.
%   CURRENTLY ONLY WORKS UP TO FLOW CONTROL ORIFICE
    for i=1:length(systable.PartName)
        if i==1
            systable.P1(i) = Ptank;
            systable.T1(i) = Ttank;
        else
            systable.P1(i) = systable.P2(i-1);
            systable.T1(i) = systable.T2(i-1);
        end
        try
            outs = valve(medium,mdot,systable.P1(i),systable.T1(i),systable.Cv(i));
            systable.P2(i) = outs(1);
            systable.T2(i) = outs(2);
        catch % If valve is choked, take choked mass flow as system mdot
            mdotout = valve(medium,mdot,systable.P1(i),systable.T1(i),systable.Cv(i));
            break
        end
    end
    if exist('mdotout','var') == 0 % If no valves are choked, orifice sets output mdot
        mdotout = chokedorifice(A,medium.gam,medium.Mw,systable.P2(end),systable.T2(end));
    end
    % If valves are choked, valve choking sets mass flow out
    delta_mdot = mdot - mdotout;
    
    if ~isempty(varargin)
        output = {delta_mdot, mdotout, systable};
    else
        output = delta_mdot;
    end
end