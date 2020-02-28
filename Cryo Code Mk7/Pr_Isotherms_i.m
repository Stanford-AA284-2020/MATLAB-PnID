% Make a pressure-density isotherm diagram.  
% This file is useful for exploring how the P_rT function behaves in various 
% regions (triple line, critical point, etc.).
% C.F. Edwards, 2-11-12 

% Provide access to support files via the Matlab path.
addpath 'Fundamental Relation Files' 
addpath 'Fundamental Relation Data'
addpath 'Setup Files' 
addpath 'Property Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2

% Put a few isotherms on the P-r plot.
Tlist = [...
    Ttrip_i(ispecies)...
    0.1*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.5*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.75*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.9*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.95*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.00*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.05*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.50*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    2.0*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.1*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    0.4*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    Tupper_i(ispecies)...
    ]

% Set limits and step size.
rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;
dr = (rmax-rmin)/steps;
% Preallocate storage...
Prisotherm  = zeros(length(Tlist),steps+1);
risotherm   = zeros(length(Tlist),steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j)
    i = 1;
    for r=rmin:dr:rmax+dr
        r;
        Prisotherm(j,i) = P_irT(ispecies,r,T);
        risotherm(j,i) = r;
        i = i+1;
    end
end

figure(1)
clf
hold on
% Put the critical point and ends of the triple line on.
plot(rcrit_i(ispecies),Pcrit_i(ispecies)/1e6,'kd')
plot([rftrip_i(ispecies) rgtrip_i(ispecies)],...
    [Ptrip_i(ispecies) Ptrip_i(ispecies)]/1e6,'ko-')
legend('Critical Point','Triple Line');
% Put the isotherms on.
for j=1:1:length(Tlist)
    plot(risotherm(j,:),Prisotherm(j,:)/1e6,'b')
end
plot(risotherm(1,:),Prisotherm(1,:)/1e6,'r')
plot(risotherm(length(Tlist),:),Prisotherm(length(Tlist),:)/1e6,'r')
hold off
xlabel('Density (kg/m^3)')
ylabel('Pressure (MPa)')
% Add some simple temperature labels.
for i=2:1:length(Tlist)-1
    text(2*rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
        num2str(Tlist(i)))
end
i = 1;
text(2*rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
    num2str(Tlist(i)),'Color','r')
i = length(Tlist);
text(2*rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
    ['T = ',num2str(Tlist(i)),' K'],'Color','r')
axis([0 rftrip_i(ispecies) -Pcrit_i(ispecies)/1e6 2*Pcrit_i(ispecies)/1e6])
% Gussy up the plot a little.
plotfixer
