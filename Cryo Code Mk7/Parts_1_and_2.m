% Make basic diagrams for O2.  
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

% Get the numerical critical point and triple line.
[Tnumcrit rnumcrit] = Critical_i(O2)
Pnumcrit = P_irT(O2,rnumcrit,Tnumcrit)
Tnumtrip = Ttrip_i(O2);
[Pnumtrip rfnumtrip rgnumtrip] = Saturation_iT(O2,Tnumtrip)
fprintf('\n')
fprintf('               T trip(K)    P trip(Pa)   rf trip(kg/m3)  rg trip(kg/m3)\n')
fprintf('-----------------------------------------------------------------------\n')
fprintf('Experimental     %5.2f        %6.2f        %6.4f        %6.6f\n',...
    Ttrip_i(O2),Ptrip_i(O2),rftrip_i(O2),rgtrip_i(O2))
fprintf('Numerical        %5.2f        %6.2f        %6.4f        %6.6f\n',...
    Tnumtrip,Pnumtrip,rfnumtrip,rgnumtrip)
fprintf('\n')
fprintf('               T crit(K)    P crit(Pa)   r crit(kg/m3)\n')
fprintf('------------------------------------------------------\n')
fprintf('Experimental    %5.4f       %5.0f      %5.5f\n',...
    Tcrit_i(O2),Pcrit_i(O2),rcrit_i(O2))
fprintf('Numerical       %5.4f       %5.0f      %5.5f\n',...
    Tnumcrit,Pnumcrit,rnumcrit)
fprintf('\n')

% Find the spinodals from the triple line up to the near-critical region.
steps = 500;
Tmax = Tnumcrit
Tmin = Ttrip_i(O2)
dT = (Tmax - Tmin)/steps;
i = 1;
% Preallocate storage:
Tplot   = zeros(steps+1);
rfsplot = zeros(steps+1);
rgsplot = zeros(steps+1);
vfsplot = zeros(steps+1);
vgsplot = zeros(steps+1);
Pfsplot = zeros(steps+1);
Pgsplot = zeros(steps+1);
rfplot  = zeros(steps+1);
rgplot  = zeros(steps+1);
vfplot  = zeros(steps+1);
vgplot  = zeros(steps+1);
Psplot  = zeros(steps+1);
musplot = zeros(steps+1);
for T=Tmin:dT:Tmax
    T;
    % Get the spinodal lines.
    rfs = Liquid_Spinodal_iT(O2,T);
    Pfs = P_irT(O2,rfs,T);
    rgs = Vapor_Spinodal_iT(O2,T);
    Pgs = P_irT(O2,rgs,T);
    
    % Save for plotting:
    Tplot(i)   = T;
    rfsplot(i) = rfs;
    rgsplot(i) = rgs;
    vfsplot(i) = 1/rfs;
    vgsplot(i) = 1/rgs;
    Pfsplot(i) = Pfs;
    Pgsplot(i) = Pgs;

    % Get the saturation lines.
    [Psat rf rg] = Saturation_iT(O2,T);
    musat = mu_irT(O2,rg,T);
    
    % Save for plotting:
    rfplot(i) = rf;
    rgplot(i) = rg;
    vfplot(i) = 1/rf;
    vgplot(i) = 1/rg;
    Psplot(i)  = Psat;
    musplot(i) = musat;

    i = i+1;
end
    
% Put a few isotherms on the P-v plot.
Tlist = [...
    Ttrip_i(O2)...
    0.1*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.25*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.5*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.75*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.9*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.95*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    1.00*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    1.05*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    1.25*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    1.50*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    2.0*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.1*(Tupper_i(O2)-Tcrit_i(O2)) + Tcrit_i(O2) ...
    0.4*(Tupper_i(O2)-Tcrit_i(O2)) + Tcrit_i(O2) ...
    Tupper_i(O2)...
    ];

% Take uniform steps in density.
rmax = rftrip_i(O2);
rmin = rgtrip_i(O2);
steps = 500;
dr = (rmax-rmin)/steps;
% Move the bounds out a bit to get the triple line.
rmin = rmin - 2*dr;
rmax = rmax + 2*dr;
dr = (rmax-rmin)/steps;
% Preallocate storage:
Prisotherm  = zeros(steps+1);
risotherm   = zeros(steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j);
    i = 1;
    for r=rmin:dr:rmax
        r;
        Prisotherm(j,i) = P_irT(O2,r,T);
        risotherm(j,i) = r;
        i = i+1;
    end
end

% Take logarithmic steps in specific volume.
steps = 500;
% Preallocate storage:
Pvlisotherm = zeros(length(Tlist),steps+1);
vlisotherm  = zeros(length(Tlist),steps+1);
Pvgisotherm = zeros(length(Tlist),steps+1);
vgisotherm  = zeros(length(Tlist),steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j)
    if(T >= Tnumcrit)
        % Split into two parts--limit at the critcal point.
        i = 1;
        % Liquid side...
        lvmax = log(1/rcrit_i(O2));
        lvmin = log(1/rftrip_i(O2));
        ldv = (lvmax-lvmin)/steps/2;
        % Nudge the bound out a bit.
        lvmin = lvmin - ldv;
        ldv = (lvmax-lvmin)/steps/2;
        for lv=lvmin:ldv:lvmax
            v = exp(lv);
            Pvlisotherm(j,i) = P_irT(O2,1/v,T);
            vlisotherm(j,i) = v;
            i = i+1;
        end
        % Vapor side...
        lvmax = log(1/rgtrip_i(O2));
        lvmin = log(1/rcrit_i(O2));
        ldv = (lvmax-lvmin)/steps/2;
        % Nudge the bound out a bit.
        lvmax = lvmax + ldv;
        ldv = (lvmax-lvmin)/steps/2;
        for lv=lvmin:ldv:lvmax
            v = exp(lv);
            Pvgisotherm(j,i) = P_irT(O2,1/v,T);
            vgisotherm(j,i) = v;
            i = i+1;
        end
    else
        % Split into two parts--limit at the spinodals.
        i = 1;
        % Liquid side...
        lvmax = log(1/Liquid_Spinodal_iT(O2,T));
        lvmin = log(1/rftrip_i(O2));
        ldv = (lvmax-lvmin)/steps/2;
        % Nudge the bound out a bit.
        lvmin = lvmin - ldv;
        ldv = (lvmax-lvmin)/steps/2;
        for lv=lvmin:ldv:lvmax
            v = exp(lv);
            Pvlisotherm(j,i) = P_irT(O2,1/v,T);
            vlisotherm(j,i) = v;
            i = i+1;
        end
        % Vapor side...
        lvmax = log(1/rgtrip_i(O2));
        lvmin = log(1/Vapor_Spinodal_iT(O2,T));
        ldv = (lvmax-lvmin)/steps/2;
        % Nudge the bound out a bit.
        lvmax = lvmax + ldv;
        ldv = (lvmax-lvmin)/steps/2;
        for lv=lvmin:ldv:lvmax
            v = exp(lv);
            Pvgisotherm(j,i) = P_irT(O2,1/v,T);
            vgisotherm(j,i) = v;
            i = i+1;
        end
    end
end

figure(1)
clf
semilogx(1/rcrit_i(O2),Pcrit_i(O2)/1e6,'bd')
hold on
semilogx(1/rnumcrit,Pnumcrit/1e6,'kd')
semilogx([1/rfnumtrip 1/rgnumtrip],[Pnumtrip/1e6 Pnumtrip/1e6],'ko')
legend('Exp. Critical','Num. Critical','Num. Triple','Location','North')
semilogx([1/rfnumtrip 1/rgnumtrip],[Pnumtrip/1e6 Pnumtrip/1e6],'k-')
for j=1:1:length(Tlist)
    semilogx(vlisotherm(j,:),Pvlisotherm(j,:)/1e6,'b')
    semilogx(vgisotherm(j,:),Pvgisotherm(j,:)/1e6,'b')
end
semilogx(vlisotherm(1,:),Pvlisotherm(1,:)/1e6,'r')
semilogx(vgisotherm(1,:),Pvgisotherm(1,:)/1e6,'r')
semilogx(vlisotherm(length(Tlist),:),Pvlisotherm(length(Tlist),:)/1e6,'r')
semilogx(vgisotherm(length(Tlist),:),Pvgisotherm(length(Tlist),:)/1e6,'r')
semilogx(vgsplot,Pgsplot/1e6,'k--',vfsplot,Pfsplot/1e6,'k--')
semilogx(vgplot,Psplot/1e6,'k',vfplot,Psplot/1e6,'k')
hold off
xlabel('Specific Volume (m^3/kg)')
ylabel('Pressure (MPa)')
for i=2:1:length(Tlist)-1
    text(2,1.4*(i)*(Pcrit_i(O2)/1e6)/(length(Tlist)),num2str(Tlist(i),'%.2f'))
end
i = 1;
text(2,1.4*(i)*(Pcrit_i(O2)/1e6)/(length(Tlist)),num2str(Tlist(i),'%.2f'),'Color','r')
i = length(Tlist);
text(2,1.4*(i)*(Pcrit_i(O2)/1e6)/(length(Tlist)),...
    ['T = ',num2str(Tlist(i),'%.0f'),' K'],'Color','r')
vmax = 1/rgtrip_i(O2);
vmin = 1/rftrip_i(O2);
scale = axis;
% axis([scale(1) scale(2) -.25 2])
plotfixer

% Make the mu-P curves.

clear Pgisotherm Plisotherm glisotherm Tlist

% Narrow the list to the vapor dome.
Tlist = [...
    Ttrip_i(O2)...
    0.1*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.25*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.5*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.75*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.9*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    0.95*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    1.00*(Tcrit_i(O2)-Ttrip_i(O2)) + Ttrip_i(O2)...
    ];

steps = 500;
% Preallocate storage:
Plisotherm = zeros(length(Tlist),steps+1);
mulisotherm = zeros(length(Tlist),steps+1);
Pgisotherm = zeros(length(Tlist),steps+1);
mugisotherm = zeros(length(Tlist),steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j)
    if(T >= Tcrit_i(O2))
        % Split into two parts--limit at the critcal point.
        i = 1;
        % Liquid side...
        rmax = rftrip_i(O2);
        rmin = rcrit_i(O2);
        dr = (rmax-rmin)/steps/2;
        % Nudge the bound out a bit.
        rmax = rmax + 2*dr;
        dr = (rmax-rmin)/steps/2;
        for r=rmin:dr:rmax
            Plisotherm(j,i)  = P_irT(O2,r,T);
            mulisotherm(j,i) = mu_irT(O2,r,T);
            i = i+1;
        end
        % Vapor side...
        rmax = rcrit_i(O2);
        rmin = rgtrip_i(O2);
        dr = (rmax-rmin)/steps/2;
        % Nudge the bound out a bit.
        rmin = rmin - 2*dr;
        dr = (rmax-rmin)/steps/2;
        for r=rmin:dr:rmax
            Pgisotherm(j,i)  = P_irT(O2,r,T);
            mugisotherm(j,i) = mu_irT(O2,r,T);
            i = i+1;
        end
    else
        % Split into two parts--limit at the spinodals.
        i = 1;
        % Liquid side...
        rmax = rftrip_i(O2);
        rmin = Liquid_Spinodal_iT(O2,T);
        dr = (rmax-rmin)/steps/2;
        % Nudge the bound out a bit.
        rmax = rmax + 2*dr;
        dr = (rmax-rmin)/steps/2;
        for r=rmin:dr:rmax
            Plisotherm(j,i)  = P_irT(O2,r,T);
            mulisotherm(j,i) = mu_irT(O2,r,T);
            i = i+1;
        end
        % Vapor side...
        rmax = Vapor_Spinodal_iT(O2,T);
        rmin = rgtrip_i(O2);
        dr = (rmax-rmin)/steps/2;
        % Nudge the bound out a bit.
        rmin = rmin - 2*dr;
        dr = (rmax-rmin)/steps/2;
        for r=rmin:dr:rmax
            Pgisotherm(j,i)  = P_irT(O2,r,T);
            mugisotherm(j,i) = mu_irT(O2,r,T);
            i = i+1;
        end
    end
end

mucrit = mu_irT(O2,rcrit_i(O2),Tcrit_i(O2))
mutrip = mu_irT(O2,rgtrip_i(O2),Ttrip_i(O2))
munumcrit = mu_irT(O2,rnumcrit,Tnumcrit)
munumtrip = mu_irT(O2,rgnumtrip,Ttrip_i(O2))

figure(2)
clf
hold on
plot(Pcrit_i(O2)/1e6,mucrit/1e6,'rd')
plot(Pnumcrit/1e6,munumcrit/1e6,'kd')
plot(Pnumtrip/1e6,munumtrip/1e6,'ko')
legend('Exp. Critical Pt.','Num. Critical Pt.','Num. Triple Pt.','Location','North')
j = 1;
plot(Plisotherm(j,:)/1e6,mulisotherm(j,:)/1e6,'r')
plot(Pgisotherm(j,:)/1e6,mugisotherm(j,:)/1e6,'r')
for j=2:1:length(Tlist)-1
    plot(Plisotherm(j,:)/1e6,mulisotherm(j,:)/1e6,'b')
    plot(Pgisotherm(j,:)/1e6,mugisotherm(j,:)/1e6,'b')
end
j = length(Tlist);
plot(Plisotherm(j,:)/1e6,mulisotherm(j,:)/1e6,'r')
plot(Pgisotherm(j,:)/1e6,mugisotherm(j,:)/1e6,'r')
plot(Psplot/1e6,musplot/1e6,'k')
hold off
xlabel('Pressure (MPa)')
ylabel('Chemical Potential (MJ/kmol)')
for i=2:1:length(Tlist)-1
    text(1.4,.5*(2-i)*(Pcrit_i(O2)/1e6)/length(Tlist),num2str(Tlist(i),'%.2f'))
end
i = 1;
text(1.4,.5*(2-i)*(Pcrit_i(O2)/1e6)/length(Tlist),...
    ['T = ',num2str(Tlist(i),'%.2f'),' K'],'Color','r')
i = length(Tlist);
text(1.4,.5*(2-i)*(Pcrit_i(O2)/1e6)/length(Tlist),num2str(Tlist(i),'%.2f'),...
    'Color','r')
scale = axis;
% axis([-.25 2 -.8 .2])
plotfixer
