function T2 = T2_interp(medium, P1, T1, P2, temp_interp)
% T2_interp Outlet Temp. Determination Through Flow Control Device
%   T2 = T2_interp(medium, P1, T1, P2, temp_interp) Calculates the
%   temperature downstream of a flow control device given the fluid medium,
%   upstream pressure and temperature, downstream pressure, and 
%   interpolation method.

if temp_interp == "isothermal"
    T2 = T1;
elseif temp_interp == "adiabatic"
    if medium.Formula == "CH4"
        p00 = -9.889E05;
        p10 = 3709;
        p01 = -0.02841;
        p20 = -1.633;
        p11 = 7.188E-05;
        p02 = -7.395E-11;
        p30 = 0.0009177;
        p21 = -4.403E-08;
        p12 = -1.345E-13;
        p03 = 7.795E-18;
        h_fit = @(P, T) p00 + (p10*T) + (p01*P) + (p20*(T.^2)) + ...
            (p11.*T.*P) + (p02*(P.^2)) + (p30*(T.^3)) + ...
            (p21.*(T.^2).*P) + (p12.*T.*(P.^2)) + (p03*(P.^3));
        h1 = h_fit(P1, T1);% Upstream Enthalpy
        h_obj = @(T2) h1 - h_fit(P2, T2);% Objective Function
        [Tll, Tul] = bracket_sign_change(h_obj, T1*0.9, T1);
        [Tll, Tul] = bisection(h_obj, Tll, Tul);
        T2 = mean([Tll, Tul]);% Result is middle of converged range
    else
        warning("Adiabatic T2 interpolation only supported for Methane, reverting to isothermal")
        T2=T1;
    end
elseif temp_interp == "isentropic"
    T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
end