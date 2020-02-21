% Provide access to support files via the Matlab path.
addpath 'Fundamental Relation Files' 
addpath 'Fundamental Relation Data'
addpath 'Setup Files' 
addpath 'Property Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n**************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;


t = 210:1:1000; p = 1:1:160;
p = p*1e5;
[T, P] = meshgrid(t, p);

sizeT = size(T); l = sizeT(1,1); b = sizeT(1,2);

r = ones(l,b); h = ones(l,b);

for i = 1:1:l
    for j = 1:1:b
        r(i,j) = rv_iTP(CH4, T(i,j), P(i,j));
        h(i,j) = h_irT(CH4, r(i,j), T(i,j));
%         r(i,j) = rv_iTP(O2, T(i,j), P(i,j));
%         h(i,j) = h_irT(O2, r(i,j), T(i,j));
    end
end

mesh(T,P/1e5,h/1000);
xlabel('Temperature [K]'); ylabel('Pressure [bar]'); zlabel('Specific Enthalpy [kJ/kg]');
plotfixer;


%%%%%%
% p00 = -8.867E05; 
% p10 = 2969;
% p01 = -0.01328;
% p20 = -0.1125;
% p11 = 1.567E-05;
% p02 = -4.197E-11;
% 
% h_fit = p00 + (p10*T) + (p01*P) + (p20*(T.^2)) + (p11.*T.*P) + (p02*(P.^2));


%%%%%% Fits from CFTool %%%%%%%%%
p00 = -9.889E05; 
p10 = 3709;
p01 = -0.02841;
p20 = -1.633;
p11 = 7.188E-05;
p02 = -7.395E-11;
p30 = 0.0009177;
p21 = -4.403E-08;
p12 = -1.345E-13;
p03 = 7.795E-18;

h_fit = p00 + (p10*T) + (p01*P) + (p20*(T.^2)) + (p11.*T.*P) + (p02*(P.^2)) + (p30*(T.^3)) + (p21.*(T.^2).*P) + (p12.*T.*(P.^2)) + (p03*(P.^3));
