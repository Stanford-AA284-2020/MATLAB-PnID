% Make contour plot of dP/dr and d2Pdr2 near critical point.  
% C.F. Edwards, 2/6/10 

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

% Get the numerical critical point data.
Tcrit_exp = Tcrit_i(N2)
rcrit_exp = rcrit_i(N2)
Pcrit_exp = Pcrit_i(N2)

% The following routine finds the numerical critical point by a numerical
% search, but this can be done by hand (graphically) as well since it is
% only one point.
[Tcrit_num rcrit_num] = Critical_i(N2)
Pcrit_num = P_irT(N2,rcrit_num,Tcrit_num)

% Make a matrix of values near the critical point.
steps = 200;
rmin = 0.9*rcrit_exp;
rmax = 1.1*rcrit_exp;
Tmin = 0.9*Tcrit_exp;
Tmax = 1.1*Tcrit_exp;
dr = (rmax-rmin)/steps;
dT = (Tmax-Tmin)/steps;
j = 1;
for T=Tmin:dT:Tmax
    Tj(j) = T;
    i = 1;
    for r=rmin:dr:rmax
        ri(i) = r;
        dPdrji(j,i) = dPdr_irT(N2,r,T);
        d2Pdr2ji(j,i) = d2Pdr2_irT(N2,r,T);
        i = i+1;
    end
    j = j+1;
end

figure(1)
clf
[C,h] = contour(ri,Tj,dPdrji/1e3,'LineWidth',2,'Color','b');
clabel(C,h)
hold on
[C,h] = contour(ri,Tj,d2Pdr2ji/1e3,'LineWidth',2,'Color',[0 .7 0]);
clabel(C,h)
plot(rcrit_exp,Tcrit_exp,'k*')
plot(rcrit_num,Tcrit_num,'r*')
plot([rmin rmax],[Tcrit_exp Tcrit_exp],'k--')
plot([rcrit_exp rcrit_exp],[Tmin Tmax],'k--')
hold off
xlabel('Density (kg/m^3)')
ylabel('Temperature (K)')
title('Oxygen')
axis([rmin rmax Tmin Tmax])
legend('dP/d\rho','d^2P/d\rho^2','CP-Exp.','CP-Num.',1)
plotfixer
