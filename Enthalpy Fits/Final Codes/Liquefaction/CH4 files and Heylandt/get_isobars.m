function isobars = get_isobars(P)
%%Returns the matrix of T-s values corresponding to the given h, and
%%pressure values listed in Plist
Setup_Props_i;
Pcrit = Pcrit_i(CH4); Ttrip = Ttrip_i(CH4);
isobars = [];
Tlist = linspace(95, 700, 20);
P

if P < Pcrit

    [Tsat rf rg] = Saturation_iP(CH4, P);
    sf = s_irT(CH4, rf, Tsat);
    sg = s_irT(CH4, rg, Tsat);

    for j = 1:length(Tlist)
        T = Tlist(j)

        if T < Tsat
            %Supercooled
            r = rl_iTP(CH4, T, P);
            s = s_irT(CH4, T, P);
            isobars = [isobars; T s];

        elseif T == Tsat
            %Saturated
            isobars = [isobars; T sf; T sg];

        else
            %Superhot
            r = rv_iTP(CH4, T, P);
            s = s_irT(CH4, T, P);
            isobars = [isobars; T s];

        end
    end
    
else

    for j = 1:length(Tlist)
        T = Tlist(j);

        %Supercooled
        r = rl_iTP(CH4, T, P);
        s = s_irT(CH4, T, P);
        isobars = [isobars; T s];

    end
end

return;
