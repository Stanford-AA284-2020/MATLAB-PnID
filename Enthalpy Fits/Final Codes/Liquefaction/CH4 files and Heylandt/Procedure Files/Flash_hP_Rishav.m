function [Tsat, y_nH2, conversion] = Flash_hp_Rishav(H, y, P_init, P_final)
% This function should return the composition of the gas mixture, and the
% percentage of liquefaction methane after constant enthalpy flash

Setup_Props_i;

if y == 0
    disp('Pure Methane');
    return;
end

[Tsat rf rg] = Saturation_iP(CH4,P_final);

r_nH2 = rv_iTP(nH2, P_final, Tsat);
h_nH2 = h_irT(nH2, r_nH2, Tsat);

% Assume total mass = 1;
M = 1;
m_nH2 = y*M;
m_CH4 = (1 - y)*M;

% Enthalpy of hydrogen
H_nH2 = m_nH2*h_nH2;

% Enthalpy of methane
H_CH4 = H - H_nH2;
h_CH4 = H_CH4/m_CH4;

% Red flags
hf_CH4 = h_irT(CH4, rf, Tsat);
hg_CH4 = h_irT(CH4, rg, Tsat);

if h_CH4 < hf_CH4 
    y_nH2 = 1; 
    conversion = 100;
    disp('Methane is in supercooled liquid state');
    return;
end

if h_CH4 > hg_CH4
    conversion = 0;
    y_nH2 = y;
    disp('Cannot flash Methane into the vapordome. Reduce enthalpy');
    return;

else

    quality = (h_CH4 - hf_CH4)/(hg_CH4 - hf_CH4);
    conversion = (1 - quality)*100;

    total_gas_mass = m_nH2 + quality*m_CH4;
    y_nH2 = m_nH2/total_gas_mass;
    
end

end


