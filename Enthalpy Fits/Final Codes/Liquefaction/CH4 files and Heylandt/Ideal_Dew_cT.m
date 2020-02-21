function [P rg rf x] = Ideal_Dew_cT(y,T,varargin)
% Return the dew-point pressure P (Pa), molefractions x, and densities (kg/m3) 
% for given vapor molefractions y and temperature T (K) as calculated by
% assuming that the liquid complement can be treated as an ideal solution.
% C.F. Edwards, 2-19-12

% No starting values are required.  The purpose of this function is to
% generate them.

global toler
global nH2 CH4 Ar
global Ru M_i

% Use information supplied if available.
switch nargin
    case 2
        % Get the inflection-point properties for the vapor composition.  
        % Doing this once will save time on rv_cTP calls.
        [Tinfl_y rinfl_y] = Pr_Inflection_c(y);
    case 4
        Tinfl_y = varargin{1};
        rinfl_y = varargin{2};
    otherwise
        disp('Incorrect number of arguments in Ideal_Dew_cT')
        return
end
Pinfl_y = P_crT(y,rinfl_y,Tinfl_y);

% Find out how many components.
N = length(y);  
if ~((N == 2)||(N == 3))
    disp('Number of components must be 2 or 3 in Ideal_Dew_cT')
    return
end

% Preallocate some storage.
x     = zeros(N,1);
xlow  = zeros(N,1);
xhigh = zeros(N,1);
Pfsi  = zeros(N,1);

% Use a small increment for P derivatives.
Pinc = 10;  % Pascals

% Find the lowest pressure at which you can make an ideal solution.
for i=1:1:length(y)
    rfs = Liquid_Spinodal_iT(i,T);
    Pfsi(i) = P_irT(i,rfs,T);
end
% The largest of these sets the lower bound.  If all are negative, use
% a small value (like 1 Pa).
Pmin = 1.01*max([Pfsi' 1]);

% Find the highest pressure at which you can evaluate the vapor mixture
% chemical potentials.
rgs = Vapor_Spinodal_cT(y,T);
Pgs = P_crT(y,rgs,T);
% Watch out for bad spinodals!
if(Pgs > Pinfl_y)
    Pmax = Pinfl_y;
else
    Pmax = Pgs;
end

% Check to see if an ideal solution is possible.  If not, return zeros.
if(Pmax < Pmin)
    disp('Pressures do not overlap in Ideal_Dew_cT')
    P = 0; rg = 0; rf = 0; x = [0 0 0];
    return
end

% Choose a starting value in this range.
Pstart = 0.5*(Pmax-Pmin) + Pmin;

% Set at starting pressure and vapor density.
P = Pstart;
rv = rv_cTP(y,T,P);
rlnH2 = rl_iTP(nH2,T,P);
rlCH4 = rl_iTP(CH4,T,P);
if(N == 3)
    % Argon is included.
    rlAr = rl_iTP(Ar,T,P);
end
    
imax = 50;
Plast = P;
for i=1:1:imax
    % Find the target chemical potentials in the vapor.
    % Note that these are the actual, mixture-model values, not an ideal
    % approximation.
    rv = rv_cTP(y,T,P,rv);
    muvCH4 = mui_icrT(CH4,y,rv,T);
    muvnH2 = mui_icrT(nH2,y,rv,T);
    
    % Find the chemical potentials of the pure liquid components.
    rlnH2 = rl_iTP(nH2,T,P,rlnH2);
    mulnH2neat = mu_irT(nH2,rlnH2,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(nH2,[1 0 0],rvnH2,T), but this is faster.
    rlCH4 = rl_iTP(CH4,T,P,rlCH4);
    mulCH4neat = mu_irT(CH4,rlCH4,T);
    % Note that you could also get the value above by using the mixture
    % function mui_icrT(CH4,[0 1 0],rvCH4,T), but this is faster.
    
    % Set the mole fractions as per an ideal solution.
    x(CH4) = exp((muvCH4 - mulCH4neat)/Ru/T);
    x(nH2) = exp((muvnH2 - mulnH2neat)/Ru/T);
    
    if(N == 3)
        % Argon is included.
        muvAr = mui_icrT(Ar,y,rv,T);
        rlAr = rl_iTP(Ar,T,P,rlAr);
        mulArneat = mu_irT(Ar,rlAr,T);
        % Note that you could also get the value above by using the mixture
        % function mui_icrT(Ar,[0 0 1],rvAr,T), but this is faster.
        x(Ar) = exp((muvAr - mulArneat)/Ru/T);
    end

    % Make the composite liquid density via Amagat.
    mnH2 = x(nH2)*M_i(nH2);
    vnH2 = mnH2/rlnH2;
    mCH4 = x(CH4)*M_i(CH4);
    vCH4 = mCH4/rlCH4;
    vmix = vnH2 + vCH4;
    mmix = mnH2 + mCH4;
    if(N == 3)
        mAr = x(Ar)*M_i(Ar);
        vAr = mAr/rlAr;
        vmix = vmix + vAr;
        mmix = mmix + mAr;
    end
    rlmix = mmix/vmix;
    
    % Check on the sum to see if it adds to unity.
    xsum = sum(x);
    f = xsum - 1;
    
    % Are we done?
    if((abs(f) < toler)&&((abs(P-Plast)/P) < toler))
        rg = rv;
        rf = rlmix;
        x = real(x/xsum);
        return
    end
    Plast = P;
    
    % Use NR to adjust pressure to get the mole fraction sum correct.
    % Choose the order of the derivative by setting the next line.
    order = 2;
    if(order == 2)
        % Use a second-order central difference for the numerical derivative.
        % Go low side.
        rv = rv_cTP(y,T,P-Pinc,rv);
        muvCH4 = mui_icrT(CH4,y,rv,T);
        muvnH2 = mui_icrT(nH2,y,rv,T);
        rl = rl_iTP(nH2,T,P-Pinc,rlnH2);
        mulnH2neat = mu_irT(nH2,rl,T);
        rl = rl_iTP(CH4,T,P-Pinc,rlCH4);
        mulCH4neat = mu_irT(CH4,rl,T);
        xlow(CH4) = exp((muvCH4 - mulCH4neat)/Ru/T);
        xlow(nH2) = exp((muvnH2 - mulnH2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P-Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xlow(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumlow = sum(xlow);
        flow = xsumlow - 1;
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvCH4 = mui_icrT(CH4,y,rv,T);
        muvnH2 = mui_icrT(nH2,y,rv,T);
        rl = rl_iTP(nH2,T,P+Pinc,rlnH2);
        mulnH2neat = mu_irT(nH2,rl,T);
        rl = rl_iTP(CH4,T,P+Pinc,rlCH4);
        mulCH4neat = mu_irT(CH4,rl,T);
        xhigh(CH4) = exp((muvCH4 - mulCH4neat)/Ru/T);
        xhigh(nH2) = exp((muvnH2 - mulnH2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P+Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xhigh(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumhigh = sum(xhigh);
        fhigh = xsumhigh - 1;
        % Form the derivative.
        dfdP = ((fhigh)-(flow))/(2*Pinc);
    else
        % Use a first-order forward difference for the numerical derivative.
        % Go high side.
        rv = rv_cTP(y,T,P+Pinc,rv);
        muvCH4 = mui_icrT(CH4,y,rv,T);
        muvnH2 = mui_icrT(nH2,y,rv,T);
        rl = rl_iTP(nH2,T,P+Pinc,rlnH2);
        mulnH2neat = mu_irT(nH2,rl,T);
        rl = rl_iTP(CH4,T,P+Pinc,rlCH4);
        mulCH4neat = mu_irT(CH4,rl,T);
        xhigh(CH4) = exp((muvCH4 - mulCH4neat)/Ru/T);
        xhigh(nH2) = exp((muvnH2 - mulnH2neat)/Ru/T);
        if(N == 3)
            % Argon is included.
            muvAr = mui_icrT(Ar,y,rv,T);
            rl = rl_iTP(Ar,T,P+Pinc,rlAr);
            mulArneat = mu_irT(Ar,rl,T);
            xhigh(Ar) = exp((muvAr - mulArneat)/Ru/T);
        end
        xsumhigh = sum(xhigh);
        fhigh = xsumhigh - 1;
        % Form the derivative.
        dfdP = ((fhigh)-(f))/(Pinc);
    end

    % Get the new pressure estimate.
    dP = -f/dfdP;
    P = P + dP;
    % Watch for out of bounds.
    if(P < Pmin)
        P = Pmin;
    end
    if(P > Pmax)
        P = Pmax;
    end
end

P = 0; rg = 0;  rf = 0; x = zeros(1,N);
disp('Fell off the end of NR loop for P in Ideal_Dew_cT')
