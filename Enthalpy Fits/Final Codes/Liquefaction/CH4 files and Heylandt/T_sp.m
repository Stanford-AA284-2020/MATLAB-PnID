function T = T_sp(y, S, P, Ti, switcher, bleed)
%%Function returns the isentropic compression state
Setup_Props_i;
y_H2 = y;
y_CH4 = 1 - y;

dT = 1; 
toler = 1.5e-3;
T = Ti;

while Ti > 0 
    
    r_nH2 = rv_iTP(nH2, T, P);
    h_nH2 = h_irT(nH2, r_nH2, T);
    r_CH4 = rv_iTP(CH4, T, P);
    h_CH4 = h_irT(CH4, r_CH4, T);
    s_nH2 = s_irT(nH2, r_nH2, T);
    s_CH4 = s_irT(CH4, r_CH4, T);

    H = y_H2*(h_nH2) + y_CH4*(h_CH4);
    S_new = y_H2*(s_nH2) + y_CH4*(s_CH4);
    S_new = S_new*bleed;
    err = abs((S_new - S)/S) 
    
    if abs((S_new - S)/S) < toler
        return;
        break;
    end
    
    if switcher == 1
        T = T + dT;
    else
        T = T - dT;
    end
    
end

end