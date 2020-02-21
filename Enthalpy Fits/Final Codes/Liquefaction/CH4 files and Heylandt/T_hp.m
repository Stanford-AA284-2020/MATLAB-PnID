function T = T_hp(y, H, P, Ti, switcher, x)
%%Function returns the isenthalpic compression state
% Switcher is used to convey the state of methane - supercooled liquid or
% superheated vapor
Setup_Props_i;
y_H2 = y;
y_CH4 = 1 - y;

dT = 0.1; 
toler = 1e-3;
T = Ti;

while Ti > 0 
    
    r_nH2 = rv_iTP(nH2, T, P);
    h_nH2 = h_irT(nH2, r_nH2, T);
    
    if switcher == 1
        r_CH4 = rv_iTP(CH4, T, P);
        h_CH4 = h_irT(CH4, r_CH4, T);
    else
        r_CH4 = rl_iTP(CH4, T, P);
        h_CH4 = h_irT(CH4, r_CH4, T);
    end
    
    s_nH2 = s_irT(nH2, r_nH2, T);
    s_CH4 = s_irT(CH4, r_CH4, T);

    H_new = y_H2*(h_nH2) + y_CH4*(h_CH4);
    H_new = H_new*x;
    err = abs((H_new - H)/H);
    
    if err < toler
        return;
        break;
    end
    
    T = T + dT;
    if T > 700
        T = NaN;
        return;
    end
end

end