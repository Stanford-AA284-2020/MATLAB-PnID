function isent = get_isenthalps(h, Plist)
%%Returns the matrix of T-s values corresponding to the given h, and
%%pressure values listed in Plist
Setup_Props_i;
Pcrit = Pcrit_i(CH4); Ttrip = Ttrip_i(CH4);
isent = [];
for i = 1:length(Plist)
    flag1 = 0; flag2 = 0; flag3 = 0;
    P = Plist(i)
    if P < Pcrit
        % Get saturation enthalpies
        [Tsat rf rg] = Saturation_iP(CH4, P);
        hf = h_irT(CH4, rf, Tsat);
        hg = h_irT(CH4, rg, Tsat);
        sf = s_irT(CH4, rf, Tsat);
        sg = s_irT(CH4, rg, Tsat);
        
        if h < hg && h > hf
            flag1 = 1
            x = (h - hf)/(hg - hf);
            s = sf + x*(sg - sf);
            T = Tsat;
            
        elseif h == hf
            s = sf;
            T = Tsat;
            
        elseif h == hg
            s = sg;
            T = Tsat;
                    
        elseif h < hf
            % Supercooled liquid
             flag2 = 1
             T = T_hp(0, h, P, Ttrip, 0, 1);
             r = rl_iTP(CH4, T, P);
             s = s_irT(CH4, r, T);
             
        else
             flag3 = 1
             T = T_hp(0, h, P, Tsat, 1, 1);
             r = rv_iTP(CH4, T, P);
             s = s_irT(CH4, r, T);
             
        end
    else
         T = T_hp(0, h, P, Ttrip, 0, 1);
         r = rl_iTP(CH4, T, P);
         s = s_irT(CH4, r, T);
    end
    
    isent = [isent; T s];
end
return;