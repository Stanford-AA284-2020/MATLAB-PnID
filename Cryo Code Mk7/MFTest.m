clc;
P1 = 70:-1:30; P1 = P1*1E5; m_dot = ones(1, length(P1));
P2 = 10*1E5;
Cd = 0.95;
A = (0.002^2)*(pi/4);

for i = 1:1:length(P1)
    T1 = 130;
    Pv = Saturation_iT(O2, T1);
    p1 = P1(1,i);
    dp = p1 - P2;
    k = (dp/(Pv - P2))^0.5;
    
    g_spi = (2*rl_iTP(O2, T1, p1)*dp)^0.5;
    
    %Find state 2
    r1 = rl_iTP(O2, T1, p1);
    s1 = s_irT(O2, r1, T1);
    h1 = h_irT(O2, r1, T1);
    s2_desired = s1; 
    T2_sat = Saturation_iP(O2, P2);
    rf2_sat = rl_iTP(O2, T2_sat, P2);
    rg2_sat = rv_iTP(O2, T2_sat, P2);
    sf2_sat = s_irT(O2, rf2_sat, T2_sat);
    sg2_sat = s_irT(O2, rg2_sat, T2_sat);
    
    if s1 < sf2_sat 
        % Remains subcooled liquid at state 2
        Thigh = T2_sat; Tlow = Ttrip_i(O2);
        sf_trip = s_irT(O2,rftrip_i(O2),Ttrip_i(O2));
        Reshigh = sf2_sat - s1;
        Reslow = sf_trip - s1;
        dT = (Thigh - Tlow)/1000;
%         T_next = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow);
%         r_next = rl_iTP(O2, T_next, P2);
%         s_next = s_irT(O2, r_next, T_next);
%         Res = s_next - s1;
        
        NFPs = 400;
        for j = 1:1:NFPs
            T_next = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow);
            r_next = rl_iTP(O2, T_next, P2);
            s_next = s_irT(O2, r_next, T_next);
            Res = s_next - s1;
            j
            if (abs(Res/s1) < toler)
                % Found the right T
                T2 = T_next;
                r2 = r_next;
                h2 = h_irT(O2, r_next, T_next);
                break;
            end
            if Res < 0
                Tlow = T_next;
                Reslow = Res;
            else
                Thigh = T_next;
                Reshigh = Res;
            end
        end
    else
        % In the vapor dome
        i
        T2 = T2_sat;
        q = (s1-sf2_sat)/(sg2_sat - sf2_sat);
        hg2_sat = h_irT(O2, rg2_sat, T2_sat);
        hf2_sat = h_irT(O2, rf2_sat, T2_sat);
        h2 = (q*(hg2_sat - hf2_sat)) + hf2_sat;
        r2 = (q*(rv_iTP(O2, T2, P2) - rl_iTP(O2, T2, P2) )) + rl_iTP(O2, T2, P2) ;
    end
    g_hem = r2*((2*(h1-h2))^0.5);
    m_dot(1,i) = Cd*A*((k*g_hem) + g_spi)/(1+k);
end

figure(1);
hold on;
plot(P1/1E5, m_dot, 'g-');
xlabel('Upstream Pressure [bar]'); ylabel('Mass Flow Rate [kg/s]');
title('Flow rate through orifice');
plotfixer;




    
    