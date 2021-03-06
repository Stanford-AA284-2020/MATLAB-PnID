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
       % Coefficients (with 95% confidence bounds):
        p00 =  -9.016e+05;%  (-9.017e+05, -9.014e+05)
        p10 =        3282;%  (3280, 3283)
        p01 =    -0.03271;%  (-0.03273, -0.03269)
        p20 =      -1.554;%  (-1.56, -1.548)
        p11 =   0.0001359;%  (0.0001358, 0.000136)
        p02 =  -2.111e-10;%  (-2.121e-10, -2.1e-10)
        p30 =     0.00258;%  (0.00257, 0.002591)
        p21 =  -2.423e-07;%  (-2.426e-07, -2.42e-07)
        p12 =   1.013e-12;%  (1.006e-12, 1.019e-12)
        p40 =  -2.057e-06;%  (-2.066e-06, -2.048e-06)
        p31 =   2.058e-10;%  (2.054e-10, 2.061e-10)
        p22 =  -1.424e-15;%  (-1.435e-15, -1.413e-15)
        p50 =   6.353e-10;%  (6.325e-10, 6.381e-10)
        p41 =  -6.706e-14;%  (-6.718e-14, -6.694e-14)
        p32 =   6.322e-19;%  (6.264e-19, 6.381e-19)
        % Linear model Poly52:
        h_fit = @(P, T) p00 + (p10*T) + (p01*P) + (p20*(T.^2)) + ...
            (p11.*T.*P) + (p02*(P.^2)) + (p30*(T.^3)) + ...
            (p21*(T.^2).*P) + (p12.*T.*(P.^2)) + (p40*T.^4) + ...
            (p31*(T.^3).*P) + (p22*(T.^2).*(P.^2)) + (p50*T.^5) + ...
            (p41*(T.^4).*P) + (p32*(T.^3).*(P.^2));
        % Goodness of fit:
        % SSE: 4.562e+09
        % R-square: 1
        % Adjusted R-square: 1
        % RMSE: 189.9
        
        h1 = h_fit(P1, T1);% Upstream Enthalpy
        h_obj = @(T2) h1 - h_fit(P2, T2);% Objective Function
        [Tll, Tul] = bracket_sign_change(h_obj, T1*0.9, T1);
        [Tll, Tul] = bisection(h_obj, Tll, Tul);
        T2 = mean([Tll, Tul]);% Result is middle of converged range
        
    elseif medium.Formula == "O2"
        % Coefficients (with 95% confidence bounds):
        p00 =  -2.368e+04;%  (-2.371e+04, -2.365e+04)
        p10 =        1144;%  (1143, 1144)
        p01 =   -0.009899;%  (-0.009902, -0.009895)
        p20 =     -0.9332;%  (-0.9343, -0.9322)
        p11 =   4.237e-05;%  (4.235e-05, 4.239e-05)
        p02 =  -3.773e-12;%  (-3.96e-12, -3.586e-12)
        p30 =    0.001746;%  (0.001744, 0.001747)
        p21 =  -7.348e-08;%  (-7.353e-08, -7.342e-08)
        p12 =   7.078e-14;%  (6.97e-14, 7.187e-14)
        p40 =   -1.34e-06;%  (-1.342e-06, -1.339e-06)
        p31 =   6.044e-11;%  (6.038e-11, 6.049e-11)
        p22 =  -1.336e-16;%  (-1.355e-16, -1.316e-16)
        p50 =   3.889e-10;%  (3.884e-10, 3.894e-10)
        p41 =  -1.916e-14;%  (-1.918e-14, -1.914e-14)
        p32 =   6.905e-20;%  (6.8e-20, 7.01e-20)
        % Linear model Poly52:
        h_fit = @(P, T) p00 + (p10*T) + (p01*P) + (p20*(T.^2)) + ...
            (p11.*T.*P) + (p02*(P.^2)) + (p30*(T.^3)) + ...
            (p21*(T.^2).*P) + (p12*T.*(P.^2)) + (p40*(T.^4)) + ...
            (p31*(T.^3).*P) + (p22*(T.^2).*(P.^2)) + (p50*(T.^5)) + ...
            (p41*(T.^4).*P) + (p32*(T.^3)*(P.^2));
        % Goodness of fit:
        % SSE: 1.463e+08
        % R-square: 1
        % Adjusted R-square: 1
        % RMSE: 34
        
        h1 = h_fit(P1, T1);% Upstream Enthalpy
        h_obj = @(T2) h1 - h_fit(P2, T2);% Objective Function
        [Tll, Tul] = bracket_sign_change(h_obj, T1*0.9, T1);
        [Tll, Tul] = bisection(h_obj, Tll, Tul);
        T2 = mean([Tll, Tul]);% Result is middle of converged range
    else
        warning("Adiabatic T2 interpolation only supported for Methane and Oxygen, reverting to isothermal")
        T2=T1;
    end
elseif temp_interp == "isentropic"
    T2 = T1*(P2/P1)^((medium.gam-1)/medium.gam);
end