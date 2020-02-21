function [Tsampled rsampled] = Vapor_Dome_c(c,Tmin,Nreal,Nsamples)
% Make a T-r array that spans the vapor dome from the vapor side to the
% liquid side.  Use Tmin as the lower temperature limit (K), 2*Nreal as the
% number of real data samples that form the basis of the interpolation, and
% Nsamples as the number of points (linear or log in density) to have
% across the dome in the final array.  The composition is fixed at value c.
% C.F. Edwards, 2-15-10

global toler
global N2 O2 Ar
global Tcrit_i rcrit_i

% The following are for ternary air only.  They are of no use with other
% compositions.
Mair = 28.9586;                 % kg/kmol
Tmaxcondentherm = 132.6312;     % K
Pmaxcondentherm = 3.78502e6;    % Pa
rmaxcondentherm = 10.4477*Mair; % kg/m3
Tmaxcondenbar   = 132.6035;     % K
Pmaxcondenbar   = 3.7891e6;     % Pa
rmaxcondenbar   = 11.0948*Mair; % kg/m3
Tcritair        = 132.5306;     % K
Pcritair        = 3.7860e6;     % Pa
rcritair        = 11.8308*Mair; % kg/m3

% See if the mixture is ternary air.  If so, we can do a little better than
% the general case.
N = length(c);
if(N == 3)
    if((c(O2) == 0.2096)&&(c(Ar) == 0.0092))
        Is_Air = 1;
        cair = c;
    else
        Is_Air = 0;
    end
else
    Is_Air = 0;
end

% Check for an essentially pure component.  If so, process separately.
pure = 0;
if(abs(1-c(N2)) < toler)
    pure = N2;
end
if(abs(1-c(O2)) < toler)
    pure = O2;
end
if(abs(1-c(Ar)) < toler)
    pure = Ar;
end

if(~pure)
    % Handle the general case.
    % Get the P-rho inflection point for scaling at the top of the dome.
    [Tinfl rinfl] = Pr_Inflection_c(c);
    Pinfl = P_crT(c,rinfl,Tinfl);
    
    % Use a small offset.
    if(Is_Air)
        Tmax = 0.99*Tinfl;
    else
        Tmax = 0.98*Tinfl;
    end

    % Use even steps up to this temperature.
    dT = (Tmax-Tmin)/(Nreal-1);

    % Get fully resolved (real) data for the two sides of the dome.
    T = Tmin;
    Tgi = zeros(1,Nreal); Pgi = zeros(1,Nreal); 
    rgi = zeros(1,Nreal); rgci = zeros(1,Nreal);
    Tfi = zeros(1,Nreal); Pfi = zeros(1,Nreal); 
    rfi = zeros(1,Nreal); rfci = zeros(1,Nreal);
    xci = zeros(N,Nreal); yci = zeros(N,Nreal);
    for i=1:1:Nreal
        T
        % Move up the vapor (dew-point) side.
        [Pdew rgdew rfdew xdew] = Dew_cT(c,T);
        % Move up the liquid (bubble-point) side.
        [Pbub rfbub rgbub ybub] = Bubble_cT(c,T);
        % Save the T-rho data.
        Tgi(i) = T;
        Pgi(i) = Pdew;
        rgi(i) = rgdew;
        xci(:,i) = xdew;
        rgci(i) = rgbub;
        Tfi(i) = T;
        Pfi(i) = Pbub;
        rfi(i) = rfbub;
        yci(:,i) = ybub;
        rfci(i) = rfdew;
        T = T+dT;
    end
else
    % Handle the pure component case.
    % Get the critical point for scaling at the top of the dome.
    Tcrit = Tcrit_i(pure);
    rcrit = rcrit_i(pure);
    Pcrit = P_irT(pure,rcrit,Tcrit);
    
    % Use a small offset.
    Tmax = 0.98*Tcrit;

    % Use even steps up to this temperature.
    dT = (Tmax-Tmin)/(Nreal-1);

    % Get fully resolved (real) data for the two sides of the dome.
    T = Tmin;
    Tgi = zeros(1,Nreal); Pgi = zeros(1,Nreal); 
    rgi = zeros(1,Nreal); rgci = zeros(1,Nreal);
    Tfi = zeros(1,Nreal+1); Pfi = zeros(1,Nreal+1); 
    rfi = zeros(1,Nreal+1); rfci = zeros(1,Nreal+1);
    xci = zeros(N,Nreal); yci = zeros(N,Nreal+1);
    for i=1:1:Nreal
        T
        % Move up the dome.
        [Psat rf rg] = Saturation_iT(pure,T);
        % Save the T-rho data.
        Tgi(i) = T;
        Pgi(i) = Psat;
        rgi(i) = rg;
        xci(:,i) = ones(N,1);
        rgci(i) = rf;
        Tfi(i) = T;
        Pfi(i) = Psat;
        rfi(i) = rf;
        yci(:,i) = ones(N,1);
        rfci(i) = rg;
        T = T+dT;
    end
    % Append the critical point to the liquid line.
    Tfi(Nreal+1)   = Tcrit;
    Pfi(Nreal+1)   = Pcrit;
    rfi(Nreal+1)   = rcrit;
    yci(:,Nreal+1) = ones(N,1);
end

% Concatenate the two to make one real data array for T and one for rho.
if(Is_Air)
    % Add the critical point between the bubble and dew data.
    Ti = [Tgi Tcritair Tfi];
    Pi = [Pgi Pcritair Pfi];
    ri = [rgi rcritair rfi];
    rci = [rfci rcritair rgci];
    cci = [xci cair yci];
else
    % We don't know the critical point in the general case.
    Ti = [Tgi Tfi];
    Pi = [Pgi Pfi];
    ri = [rgi rfi];
    rci = [rfci rgci];
    cci = [xci yci];
end

% Sample the real data.  Log or linear.  Use comments to select.
rmax = rfi(1);
rmin = rgi(1);
lrmax = log(rmax);
lrmin = log(rmin);
dlr = (lrmax-lrmin)/(Nsamples-1);
rsampled = zeros(Nsamples,1);
for i=1:1:Nsamples 
    lr = lrmin + (i-1)*dlr;
    rsampled(i) = exp(lr);
end
% Use the following line for linear spacing.
% rsampled = linspace(rgi(1),rfi(1),Nsamples);

% Find the basic spline data.
Tsampled  = interp1(ri,Ti,rsampled,'spline');
% return

% From here on is only needed if you want to make plots to check out how
% the routine is functioning.  Place a return statement before all of this,
% or comment it out, to suppress it.
Psampled  = interp1(ri,Pi,rsampled,'spline');
rcsampled = interp1(ri,rci,rsampled,'spline');
N2sampled = interp1(ri,cci(N2,:),rsampled,'spline');
O2sampled = interp1(ri,cci(O2,:),rsampled,'spline');
if(N == 3)
    Arsampled = interp1(ri,cci(Ar,:),rsampled,'spline');
end

figure(10)
clf
hold on
plot(rsampled,Tsampled,'-k')
plot(rcsampled,Tsampled,'-b')
plot(ri,Ti,'ok')
plot(rci,Ti,'ob')
if(~pure)
    plot(rinfl,Tinfl,'*r')
end
if(Is_Air)
    plot(rmaxcondentherm,Tmaxcondentherm,'dk')
    plot(rmaxcondenbar,Tmaxcondenbar,'sk')
    plot(rcritair,Tcritair,'ok')
end
hold off
xlabel('Density (kg/m^3)')
ylabel('Temperature (K)')
axis([-20 1200 60 140])
legend('Vapor Dome','Complement',1)
if(N == 3)
    title(sprintf('%.2f N2, %.2f O2, %.2f Ar',c(N2),c(O2),c(Ar)))
elseif(N==2)
    title(sprintf('%.2f N2, %.2f O2',c(N2),c(O2)))
end
plotfixer

figure(11)
clf
hold on
plot(rsampled,Psampled/1e6,'-k')
plot(rcsampled,Psampled/1e6,'-b')
plot(ri,Pi/1e6,'ok')
plot(rci,Pi/1e6,'ob')
if(~pure)
    plot(rinfl,Pinfl/1e6,'*r')
end
if(Is_Air)
    plot(rmaxcondentherm,Pmaxcondentherm/1e6,'dk')
    plot(rmaxcondenbar,Pmaxcondenbar/1e6,'sk')
    plot(rcritair,Pcritair/1e6,'ok')
end
hold off
xlabel('Density (kg/m^3)')
ylabel('Pressure (MPa)')
legend('Vapor Dome','Complement',1)
if(N == 3)
    title(sprintf('%.2f N2, %.2f O2, %.2f Ar',c(N2),c(O2),c(Ar)))
elseif(N==2)
    title(sprintf('%.2f N2, %.2f O2',c(N2),c(O2)))
end
plotfixer

figure(12)
clf
hold on
if(N == 3)
plot(rci,cci(N2,:),'o','Color',[0 .5 0])
plot(rci,cci(O2,:),'bo')
plot(rci,10*cci(Ar,:),'ro')
plot(rcsampled,N2sampled,'-','Color',[0 .5 0])
plot(rcsampled,O2sampled,'b-')
plot(rcsampled,10*Arsampled,'r-')
plot([-20 1200],[c(N2) c(N2)],'--','Color',[0 .5 0])
plot([-20 1200],[c(O2) c(O2)],'b--')
plot([-20 1200],[10*c(Ar) 10*c(Ar)],'r--')
legend('N2','O2','10*Ar',0)
else
plot(ri,cci(N2,:),'o','Color',[0 .5 0])
plot(ri,cci(O2,:),'bo')
plot(rsampled,N2sampled,'-','Color',[0 .5 0])
plot(rsampled,O2sampled,'b-')
plot([-20 1200],[c(N2) c(N2)],'--','Color',[0 .5 0])
plot([-20 1200],[c(O2) c(O2)],'b--')
legend('N2','O2',0)
end
hold off
xlabel('Complementary Phase Density (kg/m3)')
ylabel('Complementary Phase Mole Fraction')
axis([-20 1200 0 1])
if(N == 3)
    title(sprintf('%.2f N2, %.2f O2, %.2f Ar',c(N2),c(O2),c(Ar)))
elseif(N==2)
    title(sprintf('%.2f N2, %.2f O2',c(N2),c(O2)))
end
plotfixer

% Make an array to check the fit for entropy.  Since we are interested in
% making T-s or T-s-x diagrams, we need to make sure that the point spacing
% of the spline is adequater for this purpose.  (This shows why you might
% prefer log spacing to linear.)
si = zeros(1,length(ri));
for i=1:1:length(ri)
    si(i) = s_crT(c,ri(i),Ti(i));
end
ssampled = zeros(1,length(rsampled));
for i=1:1:length(rsampled)
    ssampled(i) = s_crT(c,rsampled(i),Tsampled(i));
end
if(~pure)
    sinfl = s_crT(c,rinfl,Tinfl);
end
if(Is_Air)
    smaxcondentherm = s_crT(c,rmaxcondentherm,Tmaxcondentherm);
    smaxcondenbar = s_crT(c,rmaxcondenbar,Tmaxcondenbar);
    scritair = s_crT(c,rcritair,Tcritair);
end

% Plot as a T-s diagram.
figure(13)
clf
hold on
plot(ssampled/1000,Tsampled,'k-')
plot(si/1000,Ti,'ko')
if(~pure)
    plot(sinfl/1000,Tinfl,'*r')
end
if(Is_Air)
    plot(smaxcondentherm/1000,Tmaxcondentherm,'dk')
    plot(smaxcondenbar/1000,Tmaxcondenbar,'sk')
    plot(scritair/1000,Tcritair,'ok')
end
hold off
xlabel('Entropy (kJ/kg-K)')
ylabel('Temperature (K)')
if(N == 3)
    title(sprintf('%.2f N2, %.2f O2, %.2f Ar',c(N2),c(O2),c(Ar)))
elseif(N==2)
    title(sprintf('%.2f N2, %.2f O2',c(N2),c(O2)))
end
plotfixer

% Make a P-T diagram
figure(14)
clf
hold on
plot(Tsampled,Psampled/1e6,'-k')
plot(Ti,Pi/1e6,'ok')
if(~pure)
    plot(Tinfl,Pinfl/1e6,'*r')
end
if(Is_Air)
    plot(Tmaxcondentherm,Pmaxcondentherm/1e6,'dk')
    plot(Tmaxcondenbar,Pmaxcondenbar/1e6,'sk')
    plot(Tcritair,Pcritair/1e6,'ok')
end
hold off
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
if(N == 3)
    title(sprintf('%.2f N2, %.2f O2, %.2f Ar',c(N2),c(O2),c(Ar)))
elseif(N==2)
    title(sprintf('%.2f N2, %.2f O2',c(N2),c(O2)))
end
plotfixer
