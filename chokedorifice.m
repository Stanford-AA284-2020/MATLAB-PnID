function mdot = chokedorifice(A,gam,Mw,Pt,Tt)
    % Gas Constant
    R = 8314.46261815324;% J/kmol-K
    mdot = gam/((gam+1)/2)^((gam+1)/(2*(gam-1)))*(Pt*A/sqrt(gam*R/Mw*Tt));
end