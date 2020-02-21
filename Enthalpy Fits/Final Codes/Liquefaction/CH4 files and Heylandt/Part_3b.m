% Make P-T and lnP vs. -1/T diagrams for CH4.  
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

format long
% Get the numerical critical point data.
Tcrit_exp = Tcrit_i(CH4)
rcrit_exp = rcrit_i(CH4)
Pcrit_exp = Pcrit_i(CH4)
% The following routine finds the numerical critical point by a numerical
% search, but this can be done by hand (graphically) as well since it is
% only one point.
[Tcrit_num rcrit_num] = Critical_i(CH4)
Pcrit_num = P_irT(CH4,rcrit_num,Tcrit_num)
gcrit_num = g_irT(CH4,rcrit_num,Tcrit_num)

% Get the numerical triple point data.
Ptrip_exp = Ptrip_i(CH4)
rftrip_exp = rftrip_i(CH4)
rgtrip_exp = rgtrip_i(CH4)
% The numerical triple line uses the experimental triple-point temperature
% and must satisfy saturation conditions.
[Ptrip_num rftrip_num rgtrip_num] = Saturation_iT(CH4,Ttrip_i(CH4))
gftrip_num = g_irT(CH4,rftrip_num,Ttrip_i(CH4))
ggtrip_num = g_irT(CH4,rgtrip_num,Ttrip_i(CH4))
format short

% Find the saturation line from the triple line up to the numerical
% critical point.
steps = 100;
Tmax = Tcrit_num - 5
Tmin = Ttrip_i(CH4)
dT = (Tmax - Tmin)/steps;
i = 1;
% Preallocate storage...
Tplot   = zeros(steps+1);
Psplot  = zeros(steps+1);
for T=Tmin:dT:Tmax
    T;
    % Get the saturation lines.
    [Psat rf rg] = Saturation_iT(CH4,T);
    
    % Storage for plotting:
    Tplot(i)  = T;
    Psplot(i) = Psat;
    
    i = i+1;
end

figure(1)
clf
% plot([Ttrip_i(CH4) Tcrit_i(CH4)],...
%     [Ptrip_i(CH4)/1e6 Pcrit_i(CH4)/1e6],'bo')
% hold on
% plot([Ttrip_i(CH4) Tcrit_num],...
%     [Ptrip_num/1e6 Pcrit_num/1e6],'bx')
plot(Tplot,Psplot/1e6,'b')
hold off
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
% legend('Experimental','Numerical',0)
title('P-T Plot for CH4')
plotfixer

% Set the x tick positions.
%Tticks = [13,14,15,17,20,25,30,35]'
% Set the labels.
%ticklabels = num2str(Tticks);
% Make cell array.
%Tlabels = {ticklabels};
figure(2)
clf
semilogy(-[1/Ttrip_i(CH4) 1/Tcrit_i(CH4)],...
    [Ptrip_i(CH4)/1e6 Pcrit_i(CH4)/1e6],'bo')
hold on
semilogy(-[1/Ttrip_i(CH4) 1/Tcrit_num],...
    [Ptrip_num/1e6 Pcrit_num/1e6],'bx')
semilogy(-1./Tplot,Psplot/1e6,'b')
hold off
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
legend('Experimental','Numerical',0)
scale = axis;
%axis([-1/13 -1/35 scale(3) scale(4)]) 
%set(gca,'XTick',-1 ./ Tticks)
%set(gca,'XTickLabel',Tlabels)
title('Duhring Plot for CH4')
plotfixer

