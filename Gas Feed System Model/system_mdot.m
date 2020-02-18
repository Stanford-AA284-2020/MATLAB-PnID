function [sys_mdot, systab] = system_mdot(systab, medium, mdot_test, mdot_step, mdot_coef, mdot_tol)
n_parts = length(systab{:,1});
while true
    % Propagate P & T through system for mdot_test
    mdot_out = 1;% Reset mdot_out > mdot_test to avoid choke boolean later
    for itm=1:n_parts
        % Get P1, T1 from nearest upstream element
        if itm~=1
            systab.P1(itm) = systab.P2(itm-1);
            systab.T1(itm) = systab.T2(itm-1);
        end

        % Calculate P2, T2 from P1, T1, and check for choking
        % All flow device codes assume isentropic temp relations
        if systab.Type(itm) == "valve"
            itm_out = valve(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cv(itm));
        elseif systab.Type(itm) == "regulator"
            itm_out = regulator(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cv(itm),systab.RegP2(itm),systab.RegDroop(itm));
        elseif systab.Type(itm) == "orifice"
            itm_out = orifice(medium,mdot_test,systab.P1(itm),systab.T1(itm),systab.Cd(itm),systab.A(itm)); 
        end

        % If unchoked, add P2, T2 to system table
        if length(itm_out) > 1
            systab.P2(itm) = itm_out(1);
            systab.T2(itm) = itm_out(2);
        else % if choked, get output mdot & log choke location
            mdot_out = itm_out;
            systab.Choked = repelem("N",n_parts)';
            systab.Choked(itm) = "Y";
            break
        end
    end

    % Check convergence
    if mdot_step <= mdot_tol
        sys_mdot = mdot_test; % Declare system mdot as converged
        break
    end

    % If choked, go back 1 step & refine step size to "sneak up" on choked
    % condition. If not choked, continue with same step size.
    if mdot_out < mdot_test
        mdot_test = mdot_test - mdot_step;
        mdot_step = mdot_step*mdot_coef;
    else
        mdot_test = mdot_test + mdot_step;
    end

end
