function A = orifice_size(medium,mdot,varargin)
% A = orifice_size(medium,mdot,dP,Cd,P1,T1) Calculates required
%   orifice area to produce desired pressure drop for given mass flow rate
%   and upstream fluid conditions. Express dP as fraction of upstream P.
%
% A = orifice_size(medium,mdot,P1,T1) Calculates required orifice area to 
%   produce choked flow at given mass flow rate and upstream conditions.

    if length(varargin) == 4
        dP = varargin{1};
        Cd = varargin{2};
        P1 = varargin{3};
        T1 = varargin{4};

        rho1 = PREoS(medium,"rho",P1,T1);
        A = mdot/(Cd*sqrt(2*rho1*(dP*P1)));
    else
        P1 = varargin{1};
        T1 = varargin{2};
        k = medium.gam;
        Mw = medium.Mw;
        R = 8314.46261815324;% J/kmol-K, Gas Constant
        
        k_coeff = k/((k+1)/2)^((k+1)/(2*(k-1)));
        % mdot = k_coeff*P1*A/sqrt(k*R/Mw*T1);
        A = mdot*sqrt(k*R/Mw*T1)/(P1*k_coeff);
    end
end