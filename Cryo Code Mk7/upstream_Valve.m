function [P_up T_up] = upstream_Valve(ispecies, P_down,T_down,Cd,Area,mdot,Pmax)
% Return the pressure and temperature upstream of a valve for the given
% downstream conditions, mass flow rate and the outlet area

global Tcrit_i rcrit_i Ttrip_i Tupper_i
global toler
persistent inumcrit Tnumcrit rnumcrit Pnumcrit

% The first time through you must set the species.
if(isempty(inumcrit))
    inumcrit = ispecies;
end

% Get numerical critical point if not already available.
if(isempty(Pnumcrit)||(inumcrit ~= ispecies))
    inumcrit = ispecies;
    try
        % This will work if you have Critical_i
        [Tnumcrit rnumcrit] = Critical_i(ispecies);
    catch
        % If you don't, assume its oxygen
        disp('Species must be ispecies in rT_ihP')
        Tnumcrit = 154.5994;
        rnumcrit = 426.9340;
     end
    Pnumcrit = P_irT(inumcrit,rnumcrit,Tnumcrit);
end

% Shorten the names.
Tcrit  = Tcrit_i(ispecies);
rcrit  = rcrit_i(ispecies);
Ttrip  = Ttrip_i(ispecies);
Tupper = Tupper_i(ispecies);
hcrit  = h_irT(ispecies,rcrit,Tcrit);


% Calculate enthalpy for the downstream state
T_sat_down = Saturation_iP(ispecies, P_down);

if T_down < T_sat_down
    % Subcooled Liquid
    r_down = rl_iTP(ispecies, T_down, P_down);
    h_down = h_irT(ispecies, r_down, T_down);
    
elseif T_down > T_sat_down
    % Superheated Vapor
    r_down = rv_iTP(ispecies, T_down, P_down);
    h_down = h_irT(ispecies, r_down, T_down);
    
else
    % Inside the vapordome
    rl_down = rl_iTP(ispecies, T_down, P_down);
    rv_down = rv_iTP(ispecies, T_down, P_down);
    hl_down = h_irT(ispecies, rl_down, T_down);
    hv_down = h_irT(ispecies, rv_down, T_down);
    q = 0;
    v = 1/rl_down + q*(1/rv_down - 1/rl_down);
    r_down = 1/v;
    h_down = hl_down + q*(hv_down - hl_down);
end

    
LHS = (mdot/(Area*Cd))^2;
Plow = P_down*1.01; % To avoid unreasonable answers
Phigh = Pmax*2; %Tank initial pressure

[r_high, T_high] = rT_ihP(ispecies, h_down, Phigh);
[r_low, T_low] = rT_ihP(ispecies, h_down, Plow);

RHS_high = 2*r_high*(Phigh - P_down);
RHS_low = 2*r_low*(Plow - P_down);

Reshigh = RHS_high - LHS;
Reslow = RHS_low - LHS;


NFPs = 200;
for j=1:1:NFPs

    P_next = Plow + (Phigh - Plow)*(-Reslow)/(Reshigh - Reslow);
    [r_next, T_next] = rT_ihP(ispecies, h_down, P_next);
    RHS_next = 2*r_next*(P_next - P_down);
    Res   = RHS_next - LHS;

    % Check for convergence. Total Volume should be V0
    if(abs(Res/LHS) < toler) 
        P_up = P_next;          % Found the solution to tolerance.
        T_up = T_next;
        break;
    end
    if(Res > 0) % Assumed Pressure is higher
        Phigh   = P_next;    % Midpoint data is on the high side of P.
        Reshigh = Res;
    else
        Plow    = P_next;    % Midpoint data is on the low side of P.
        Reslow  = Res;
    end
end
end

