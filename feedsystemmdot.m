function [delta_mdot, orificemdot, systable] = feedsystemmdot(systable, medium, Ptank, Ttank, A, mdot)
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
        [systable.P2(i), systable.T2(i)] = valve(medium,mdot,systable.P1(i),systable.T1(i),systable.Cv(i));
    end
    orificemdot = chokedorifice(A,medium.gam,medium.Mw,systable.P2(end),systable.T2(end));
    delta_mdot = mdot - orificemdot;
end