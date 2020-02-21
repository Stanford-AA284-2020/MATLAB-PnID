function [T s] = Liquid_Isobar_cPTN(c,P,Tmin,Tmax,Nsamples)
% Return T (K) and s (J/kg-K) arrays at pressure P (Pa) that span from Tmin
% to Tmax (K)

dT = (Tmax-Tmin)/(Nsamples-1);
T = zeros(Nsamples,1); s = zeros(Nsamples,1);
for i=1:1:Nsamples
    T(i) = Tmin + (i-1)*dT;
    r = rl_cTP(c,T(i),P);
    s(i) = s_crT(c,r,T(i));
end