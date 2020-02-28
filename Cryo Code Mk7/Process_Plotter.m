 load results.mat;
 Setup_Props_i;
 ispecies = O2;
 index = 51;
 
 
global Tcrit_i rcrit_i Ttrip_i Tupper_i
global toler
% persistent inumcrit Tnumcrit rnumcrit Pnumcrit

%% Defining the initial conditions %%

% Get the critical and triple data for Oxygen
Tc_O2 = Tcrit_i(O2);
rc_O2 = rcrit_i(O2);
Pc_O2 = Pcrit_i(O2);
sc_O2 = s_irT(O2,rc_O2,Tc_O2);
Tt_O2 = Ttrip_i(O2);
Pt_O2 = Ptrip_i(O2);
sft_O2 = s_irT(O2,rftrip_i(O2),Tt_O2);
sgt_O2 = s_irT(O2,rgtrip_i(O2),Tt_O2);
 
 P1 = P(1, index);
 P_global = downstream_P(index, :);
 T_global = downstream_T(index, :);
 P2 = P_global(1,1); P3 = P_global(1,2); P4 = P_global(1,3); P5 = P_global(1,4);
 T2 = T_global(1,1); T3 = T_global(1,2); T4 = T_global(1,3); T5 = T_global(1,4);
 
 T1 = T_O2(1, index);
 r1 = rl_iTP(ispecies, T1, P1);
 r2 = rl_iTP(ispecies, T2, P2);
 r3 = rl_iTP(ispecies, T3, P3);
 r4 = rl_iTP(ispecies, T4, P4);
 r5 = rl_iTP(ispecies, T5, P5);
 Pthroat = 0.9E6;
%Find state at the venturi throat
s1 = s_irT(ispecies, r1, T1);
h1 = h_irT(ispecies, r1, T1);
sthroat_desired = s1; 
Tthroat_sat = Saturation_iP(ispecies, Pthroat);
rfthroat_sat = rl_iTP(ispecies, Tthroat_sat, Pthroat);
rgthroat_sat = rv_iTP(ispecies, Tthroat_sat, Pthroat);
sfthroat_sat = s_irT(ispecies, rfthroat_sat, Tthroat_sat);
sgthroat_sat = s_irT(ispecies, rgthroat_sat, Tthroat_sat);

if s1 < sfthroat_sat 
    % Remains subcooled liquid at state throat
    Thigh = Tthroat_sat; Tlow = Ttrip_i(ispecies);
    sf_trip = s_irT(ispecies,rftrip_i(ispecies),Ttrip_i(ispecies));
    Reshigh = sfthroat_sat - s1;
    Reslow = sf_trip - s1;
    dT = (Thigh - Tlow)/1000;

    NFPs = 400;
    for j = 1:1:NFPs
        T_next = Tlow + (Thigh - Tlow)*(-Reslow)/(Reshigh - Reslow);
        r_next = rl_iTP(ispecies, T_next, Pthroat);
        s_next = s_irT(ispecies, r_next, T_next);
        Res = s_next - s1;
        j
        if (abs(Res/s1) < toler)
            % Found the right T
            Tthroat = T_next;
            rthroat = r_next;
            hthroat = h_irT(ispecies, r_next, T_next);
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
    Tthroat = Tthroat_sat;
    q = (s1-sfthroat_sat)/(sgthroat_sat - sfthroat_sat);
    hgthroat_sat = h_irT(ispecies, rgthroat_sat, Tthroat_sat);
    hfthroat_sat = h_irT(ispecies, rfthroat_sat, Tthroat_sat);
    hthroat = (q*(hgthroat_sat - hfthroat_sat)) + hfthroat_sat;
    sthroat = (q*(sgthroat_sat - sfthroat_sat)) + sfthroat_sat;
    rthroat = (q*(rv_iTP(ispecies, Tthroat, Pthroat) - rl_iTP(ispecies, Tthroat, Pthroat) )) + rl_iTP(ispecies, Tthroat, Pthroat) ;
end

% Find entropy, enthalpy, specific volume at states 2, 3, 4, 5
s2 = s_irT(ispecies, r2, T2); h2 = h_irT(ispecies, r2, T2);
s3 = s_irT(ispecies, r3, T3); h3 = h_irT(ispecies, r3, T3);
s4 = s_irT(ispecies, r4, T4); h4 = h_irT(ispecies, r4, T4);

T5_sat = Saturation_iP(ispecies, P5);
if T5 < T5_sat
    r5 = rl_iTP(ispecies, T5, P5);
    s5 = s_irT(ispecies, r5, T5); 
    h5 = h_irT(ispecies, r5, T5);
elseif T5 > T5_sat
    r5 = rv_iTP(ispecies, T5, P5);
    s5 = s_irT(ispecies, r5, T5); 
    h5 = h_irT(ispecies, r5, T5);
else
    rf5 = rl_iTP(ispecies, T5, P5);
    rg5 = rv_iTP(ispecies, T5, P5);
    sf5 = s_irT(ispecies, rf5, T5); 
    hf5 = h_irT(ispecies, rf5, T5);
    sg5 = s_irT(ispecies, rg5, T5); 
    hg5 = h_irT(ispecies, rg5, T5);
    q = (h4 - hf5)/(hg5 - hf5);
    h5 = h4;
    r5 = rf5 + q*(rg5 - rf5);
    s5 = sf5 + q*(sg5 - sf5);
end

T_Vector = [T1; Tthroat; T2; T3; T4; T5];
P_Vector = [P1; Pthroat; P2; P3; P4; P5];
v_Vector = 1./[r1; rthroat; r2; r3; r4; r5];
h_Vector = [h1; hthroat; h2; h3; h4; h5];
s_Vector = [s1; sthroat; s2; s3; s4; s5];

figure(1);
hold on;
semilogx(v_Vector, P_Vector./1E5, 'b--o');
plotfixer;

figure(2);
clf;
plot(s_Vector, T_Vector, 'bo');
plotfixer;

%% Coordinates for PV
P_coordinates = P_Vector;
P_coordinates(3,1) = P_coordinates(3,1)*1.06;
P_coordinates(4,1) = P_coordinates(4,1)*1.00;
P_coordinates(5,1) = P_coordinates(5,1)*0.94;
V_coordinates = v_Vector*1.05;
state_labels = {'1. Tank', '2. Venturi Throat', '3. CK', '4. BV', '5. INJ', '6. CH'};
text(V_coordinates, P_coordinates./1E5, state_labels, 'Color', [0, 0, 1]);
plotfixer;

 Tlist = [90 100 110 120 130 140 180 200 250 300 500 800 1000];
 for i = 1:1:length(Tlist)
     text_T{i} = [num2str(Tlist(i)), ' K'];
 end     
 
 P_IsoTcoordinates = [60 50 7 12 19 29 60 60 60 60 60 65 70];
 V_IsoTcoordinates = [8E-4 7.8E-4 6E-3 6E-3 6E-3 6E-3 6E-3 8E-3 1.1E-2 1.5E-2 2.5E-2 3.5E-2 4.2E-2];
 text(V_IsoTcoordinates, P_IsoTcoordinates, text_T, 'Color', [1, 0, 0]);
 plotfixer;