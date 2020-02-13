function out = PREoS(medium, find, in1, in2)
    R = 8314.46261815324;% J/kmol-K, Gas Constant
    if medium == "Methane"
        % Formula: CH4
        % CAS No.: 74828
        Mw = 16.043;% g/mol
        Tc = 190.564;% K
        Pc = 4.59e6;% Pa
        Vc = 0.099;% m3/kmol
        Zc = 0.286;
        w = 0.011;% Acentric Factor
    elseif medium == "Oxygen"
        % Formula: O2
        % CAS No.: 7782447
        Mw = 31.999;% kg/kmol
        Tc = 154.58;% K
        Pc = 5.02e6;% Pa
        Vc = 0.074;% m3/kmol
        Zc = 0.287;
        w = 0.020;% Acentric Factor
    elseif medium == "Nitrogen"
        % Formula: N2
        % CAS No.: 7727379
        Mw = 28.014;% kg/kmol
        Tc = 126.2;% K
        Pc = 3.39e6;% Pa
        Vc = 0.089;% m3/kmol
        Zc = 0.288;
        w = 0.037;% Acentric Factor
    elseif medium == "Air"
        % Formula: 
        % CAS No.: 132259100
        Mw = 28.951;% kg/kmol
        Tc = 132.45;% K
        Pc = 3.79e6;% Pa
        Vc = 0.092;% m3/kmol
        Zc = 0.318;
        w = 0.000;% Acentric Factor
    end
    a_Tc = 0.45724*R^2*Tc^2/Pc;% a at critical temp
    b = 0.07780*R*Tc/Pc;% b at critical temp = b at desired temp
    kap = 0.37464 + 1.54226*w - 0.26992*w^2;% Characteristic constant
    
    % Interpret equation form
    if find == "P"
        T = in1;% K
        rho = in2;% kg/m^3
        Vm = Mw/rho;% m^3/kmol
        % Use typical form (Eq 4) of PREoS to find P explicitly
        P = (R*T)/(Vm - b) - a(T)/(Vm*(Vm+b) + b*(Vm - b));
        out = P;
    elseif find == "rho"
        P = in1;% Pa
        T = in2;% K
        A = a(T)*P/(R^2*T^2);
        B = b*P/(R*T);
        % Use cubic form of PREoS to find compressibility factor
        Z = roots([1; B-1; A-3*B^2-2*B; B^3+B^2-A*B]);
        Z = Z(imag(Z)==0); % Discard imaginary roots
        out = Z;
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
        preos = @(T) (R*T)/(Vm - b) - a(T)/(Vm*(Vm+b) + b*(Vm - b)) - P;
        T = fzero(preos,Tc);
        out = T;
    end
    
    function a_T = a(T)
        alpha = (1 + kap*(1-sqrt(T/Tc)))^2;% Dimensionless fxn of reduced temp
        a_T = a_Tc*alpha;% a at desired temp
    end
end