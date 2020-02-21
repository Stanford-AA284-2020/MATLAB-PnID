% Loads data into an array for dealing with multiple species.
% Data file for methane properties using the data from 
% "Thermodynamic Properties of Cyrogenic Fluids" by R.T. Jacobsen,
% S.G. Penoncello, and E.W. Lemmon, Plenum Press, 1997.
% Code built off of CFE 1/23/08 version for CH4 
% Mohammed Ayub and Rishav Choudhary 5/19/2017

%Concerns: 
%1: Coefficient of FR coefficients and exponents has missing terms. These
%have been replaced with NaN, though I am not sure if that's the best
%approach - Ayub----- No, its not! ---- Rishav

% Constants and ideal-gas ref. state properties: Updated ------ 05/21/17
Ru          = 8314.34;                          % J/kmol-K
M_i(CH4)    = 16.0428;                          % kg/kmol
R_i(CH4)    = Ru/M_i(CH4);                      % J/kg-K
wP_i(CH4)   = 0.0086;                           % Pitzer's acentric factor
Tref_i(CH4) = 298.15;                           % K
Pref_i(CH4) = 0.101325e6;                          % Pa
rref_i(CH4) = Pref_i(CH4)/R_i(CH4)/Tref_i(CH4); % kg/m3
href = 0/M_i(CH4);                        % J/kg
sref = 0/M_i(CH4);                       % J/kg-K



% Fixed-point properties and limits: 
Tcrit_i(CH4) = 190.564;               % K
Pcrit_i(CH4) = 4.5992e6;            % Pa
rcrit_i(CH4) = 162.66;     % kg/m3

delta0 = rref_i(CH4)/rcrit_i(CH4);
tau0 = 1/(Tref_i(CH4)/Tcrit_i(CH4));

Ttrip_i(CH4)  = 90.6941;              % K
Ptrip_i(CH4)  = 0.011696e6;         % Pa
rftrip_i(CH4) = 451.48;             % kg/m3 -----Needs Update --- Updated @ 5 PM, 05/21/17
rgtrip_i(CH4) = 0.25074;            % kg/m3 ------ Needs Update --- Updated @ 5 PM, 05/21/17

Tupper_i(CH4) = 625;                % K
Tlower_i(CH4) = 90.6941;              % K
Pupper_i(CH4) = 1000e6;               % Pa
rupper_i(CH4) = 36.2029*M_i(CH4);   % kg/m3

% Ideal Cp_s - Refer to table 5-58
% There are 5 einstein terms and a constant term
N1 = 5.46;
N1001 = 0.8449e-2;
N1002 = 4.6942;
N1003 = 3.4865;
N1004 = 1.6572;
N1005 = 1.4115;

N_k = [N1001; N1002; N1003; N1004; N1005];

A1001 = 648;
A1002 = 1957;
A1003 = 3895;
A1004 = 5705;
A1005 = 15080;

A_k = [A1001; A1002; A1003; A1004; A1005];

% Translating these into the format of the a0 function gives: 
FR_Npoly0(CH4) = 2;
FR_Neinst(CH4) = 5; % Should be 5 ----- error in the input values for calculation of Einstein terms
FR_Nspec0(CH4) = 0;
FR_t0(1,CH4)    = 0;
FR_t0(2,CH4)    = 1;
FR_N0(5,CH4)  = N1001/Tcrit_i(CH4);
FR_N0(6,CH4)  = N1002/Tcrit_i(CH4);
FR_N0(7,CH4)  = N1003/Tcrit_i(CH4);
FR_N0(8,CH4)  = N1004/Tcrit_i(CH4);
FR_N0(9,CH4)  = N1005/Tcrit_i(CH4);
FR_N0(3,CH4)   = N1-1;
FR_N0(4,CH4)   = 0;
FR_gamma0(5,CH4) = A1001/Tcrit_i(CH4);
FR_gamma0(6,CH4) = A1002/Tcrit_i(CH4);
FR_gamma0(7,CH4) = A1003/Tcrit_i(CH4);
FR_gamma0(8,CH4) = A1004/Tcrit_i(CH4);
FR_gamma0(9,CH4) = A1005/Tcrit_i(CH4);

% The constants associated with the zero and first power terms in tau
% set the reference state energy and entropy.  These must be adjusted to
% give the desired state.
% FR_N0(2,CH4) = 0;    % Use this to set h0ref first.
% FR_N0(1,CH4) = 0;    % Then use this to set s0ref (depends on above).
FR_N0_1_term1 = - 1 + log(tau0/delta0);
FR_N0_1_term2 = N1*(1 - log(tau0));
FR_N0_1_term3 = 0;
FR_N0_1_term4 = 0;
FR_N0_1_term5 = 0;

for i = 1:5
    FR_N0_1_term3 = FR_N0_1_term3 + (N_k(i,1)*A_k(i,1)/(exp(A_k(i,1)*tau0/Tcrit_i(CH4)) - 1))*(tau0/Tcrit_i(CH4));
    FR_N0_1_term4 = FR_N0_1_term4 + (N_k(i,1)*A_k(i,1))*(tau0/Tcrit_i(CH4));
    FR_N0_1_term5 = FR_N0_1_term5 - (N_k(i,1)*(log(-1 + exp(A_k(i,1)*tau0/Tcrit_i(CH4)))));
end

FR_N0_2_term1 =-N1/tau0;
FR_N0_2_term2 = 0;
FR_N0_2_term3 = 0;

for i = 1:5
    FR_N0_2_term2 = FR_N0_2_term2 - (N_k(i,1)*A_k(i,1)/(exp(A_k(i,1)*tau0/Tcrit_i(CH4)) - 1))*(1/Tcrit_i(CH4));
%     FR_N0_2_term3 = FR_N0_2_term3 - (N_k(i,1)*A_k(i,1))*(1/Tcrit_i(CH4));
end

FR_N0(2,CH4) =  FR_N0_2_term1 + FR_N0_2_term2 + FR_N0_2_term3;    % h0ref is 0
FR_N0(1,CH4) =  FR_N0_1_term1 + FR_N0_1_term2 + FR_N0_1_term3 + FR_N0_1_term4 + FR_N0_1_term5;   % Then use this to set s0ref (depends on above) -------------- Needs update
%FR_N0(1,CH4) = 0;
%FR_N0(2,CH4) = 0;



% Number of terms of each type in the residual part of the FR:
FR_Npoly(CH4)  = 22; 
FR_Nexp(CH4)   = 0;
FR_Ngaus(CH4)  = 40-22;



% Exponents & Coefficients from table 5.59 of the book
Table_5_59 = [
% 1	0.4367901028e-1     1	-0.5    0	0
% 2	0.6709236199        1	0.5     0	0
% 3	-1.765577859        1	1       0	0
% 4	0.8582330241        2	0.5     0	0
% 5	-1.206513052        2	1       0	0
% 6	0.512046722         2	1.5     0	0
% 7	-0.4000010791e-3    2	4.5     0	0
% 8	-0.1247842423e-1	3	0       0	0
% 9	0.3100269701e-1     4	1       0	0
% 10	0.1754748522e-2     4	3       0	0
% 11	-0.3171921605e-5	8	1       0	0
% 12	-0.224034684e-5     9	3       0	0
% 13	0.2947056156e-6     10	3       0	0
% 14	0.1830487909        1	0       1	1
% 15	0.1511883679        1	1       1	1
% 16	-0.4289363877       1	2       1	1
% 17	0.6894002446e-1     2	0       1	1
% 18	-0.1408313996e-1	4	0       1	1
% 19	-0.306305483e-1     5	2       1	1
% 20	-0.2969906708e-1	6	2       1	1
% 21	-0.1932040831e-1	1	5       2	1
% 22	-0.1105739959       2	5       2	1
% 23	0.9952548995e-1     3	5       2	1
% 24	0.8548437825e-2     4	2       2	1
% 25	-0.6150555662e-1	4	4       2	1
% 26	-0.4291792423e-1	3	12      3	1
% 27	-0.181320729e-1     5	8       3	1
% 28	0.344590476e-1      5	10      3	1
% 29	-0.238591945e-2     8	10      3	1
% 30	-0.1159094939e-1	2	10      4	1
% 31	0.6641693602e-3     3	14      4	1
% 32	-0.237154959e-1     4	12      4	1
% 33	-0.3961624905e-1	4	18      4	1
% 34	-0.1387292044e-1	4	22      4	1
% 35	0.3389489599e-1     5	18      4	1
% 36	-0.2927378753e-2	6	14      4	1
% 37	0.9324799946e-4     2	2       0	0
% 38	-6.287171518        0	0       0	0
% 39	0.1271069467e2      0	1       0	0
% 40	-6.423953466        0	2       0	0
%  ];
1    -8.19267257565         0   3   0   0
2     0.130546404892        0   4   0   0
3    -0.11540517276e-1      0   5   0   0
4     0.839918647516e-1     1   0   0   0
5     1.33747820573         1  0.5  0   0
6    -2.79507010112         1   1   0   0
7     0.482999102128        1   2   0   0
8    -0.148739103766        1   3   0   0
9     0.25598353866e-1      2   0   0   0
10   -0.541546874713e-1     2   1   0   0
11    0.37223474567         2   2   0   0
12    1.8938236094          2   3   0   0
13    0.267165188522e-2     3   0   0   0
14    0.116573604971        3   1   0   0
15   -0.273795055848        3   2   0   0
16   -0.103129506212e-1     4   1   0   0
17    0.397065461542e-1     5   2   0   0
18   -0.840274234789e-1     5   3   0   0
19   -0.110411368495e-1     6   2   0   0
20    0.104800406095e-2     7   2   0   0
21    0.902088706725e-2     7   3   0   0
22   -0.129256450195e-2     8   3   0   0
23    8.19267257565         0   3   2   0.91464479
24   -0.130546404892        0   4   2   0.91464479
25    0.11540517276e-1      0   5   2   0.91464479
26    5.66493981456         2   3   2   0.91464479
27   -0.297031070059        2   4   2   0.91464479
28    0.105554740466e-1     2   5   2   0.91464479
29    1.80209046818         4   3   2   0.91464479
30   -0.135838960943        4   4   2   0.91464479
31    0.394549979521e-1     4   5   2   0.91464479
32    0.353263983498        6   3   2   0.91464479
33   -0.192814451572e-1     6   4   2   0.91464479
34    0.120291028247e-1     6   5   2   0.91464479
35    0.416535430489e-1     8   3   2   0.91464479
36   -0.440891835846e-2     8   4   2   0.91464479
37   -0.980954445376e-3     8   5   2   0.91464479
38    0.53460864973e-2      10  3   2   0.91464479
39    0.275654386288e-3     10  4   2   0.91464479
40   -0.179444975323e-3     10  5   2   0.91464479
 ];
% Note the we are using the Gaussian terms to implement the gamma-weighted
% exponential terms.  As such, we do not need the 5th column above, and the
% gamma coefficient becomes eta below (with the other parts of the Gaussian
% term set to zero).
FR_N(1:40,CH4)       = Table_5_59(:,2);
FR_d(1:40,CH4)       = Table_5_59(:,3);
FR_t(1:40,CH4)       = Table_5_59(:,4);
FR_eta(1:40,CH4)     = Table_5_59(:,6);
FR_beta(1:40,CH4)    = zeros(40,1);
FR_gamma(1:40,CH4)   = zeros(40,1);
FR_epsilon(1:40,CH4) = zeros(40,1);

% % Updates --- extra Gaussian terms for Methane
% % Update Betas first
% FR_beta(37,CH4) =  200;
% FR_beta(38,CH4) =  250;
% FR_beta(39,CH4) =  250;
% FR_beta(40,CH4) =  250;
% 
% % Update Gammas now (CFE uses gamma instead of delta - the term used in the
% % book)
% FR_gamma(37,CH4) = 1.07;
% FR_gamma(38,CH4) = 1.11;
% FR_gamma(39,CH4) = 1.11;
% FR_gamma(40,CH4) = 1.11;
% 
% % Update Epsilions now (all 1)
% FR_epsilon(37,CH4) = 1;
% FR_epsilon(38,CH4) = 1;
% FR_epsilon(39,CH4) = 1;
% FR_epsilon(40,CH4) = 1;
% 
% % Finally, Update Etas
% FR_eta(37,CH4) = 20;
% FR_eta(38,CH4) = 40;
% FR_eta(39,CH4) = 40;
% FR_eta(40,CH4) = 40;



