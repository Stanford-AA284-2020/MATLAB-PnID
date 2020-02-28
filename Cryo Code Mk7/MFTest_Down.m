clc; close all;
P2 = 40:-0.5:1; P2 = P2*1E5; m_dot = ones(1,length(P2)); G_HEM = ones(1,length(P2)); R2 = ones(1,length(P2));
G_SPI = ones(1,length(P2));
P1 = 45*1E5;
Cd = 0.66;
A = (0.002^2)*(pi/4);

for i = 1:1:length(P2)
    T1 = 132;
%     Pv = Saturation_iT(O2, T1);
    p2 = P2(1,i);
    dp = P1 - p2;
%     k = (dp/(Pv - p2))^0.5;
    
    g_spi = (2*rl_iTP(O2, T1, P1)*dp)^0.5;
    
    %Find state 2
    r1 = rl_iTP(O2, T1, P1);
    s1 = s_irT(O2, r1, T1);
    h1 = h_irT(O2, r1, T1);
    s2_desired = s1; 
    T2_sat = Saturation_iP(O2, p2);
    rf2_sat = rl_iTP(O2, T2_sat, p2);
    rg2_sat = rv_iTP(O2, T2_sat, p2);
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
            T_next = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow)
            r_next = rl_iTP(O2, T_next, p2);
            s_next = s_irT(O2, r_next, T_next);
            Res = s_next - s1;
            
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
        q = (s1-sf2_sat)/(sg2_sat - sf2_sat)
        hg2_sat = h_irT(O2, rg2_sat, T2_sat);
        hf2_sat = h_irT(O2, rf2_sat, T2_sat);
        h2 = (q*(hg2_sat - hf2_sat)) + hf2_sat;
        r2 = (q*(rv_iTP(O2, T2, p2) - rl_iTP(O2, T2, p2) )) + rl_iTP(O2, T2, p2) ;
%         r2 = (1-q)*rl_iTP(O2, T2, p2) ;
%         h2 = (1-q)*hf2_sat ;
        Pv = Saturation_iT(O2, T1);
        k = (dp/(Pv - p2))^0.5;
    end
    g_hem = r2*((2*(h1-h2))^0.5); 
    m_dot(1,i) = Cd*A*((k*g_hem) + g_spi)/(1+k);
    G_HEM(1,i) = g_hem; R2(1,i) = r2; G_SPI(1,i) = g_spi;
end

figure(1);
hold on;
plot(P2/1E5, m_dot, 'k-');
xlabel('Downstream Pressure [bar]'); ylabel('Mass Flow Rate [kg/s]');
title('Flow rate through orifice');
plotfixer;

figure(2);
hold on;
plot(P2/1E5, G_HEM, 'r-');
plot(P2/1E5, G_SPI, 'k-');
xlabel('Downstream Pressure [bar]'); ylabel('Mass Flux HEM [kg/m^2-s]');
title('Homogenous Equilibrium Mass Flux');
plotfixer;

figure(3);
hold on;
plot(P2/1E5, R2, 'b-');
xlabel('Downstream Pressure [bar]'); ylabel('Density [kg/m^3]');
title('Downstream Density');
plotfixer;




    
    