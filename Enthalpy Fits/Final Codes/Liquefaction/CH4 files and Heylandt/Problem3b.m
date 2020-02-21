%% Problem 3b

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

clear all
close all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N     = 2;
Setup_Air_Props;

% Set fixed-point values that are specific to air.  See Lemmon et al. 2000.
% The composition is: N2:O2:Ar = 0.7812:0.2096:0.0092
% Be sure to change these as needed if you use engineering air (79:21).
% And ignore these if you use an arbitrary composition.  They are for air.
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
% Bottom of dome conditions for air:
Tsolidair       = 59.75;        % K
Psolidair       = 5265;         % Pa
rsolidair       = 33.067*Mair;  % kg/m3
% Upper limit conditions for the model:
Tupper          = 870;          % K
Pupper          = 100e6;        % Pa
rupper          = rsolidair;
% Lower limit to stay just above solid air:
Tlower          = 60;           % K

% Set an air composition.
c = zeros(N,1);
if(N == 3)
    % Use ternary air.
    c(O2) = 0.2096;
    c(Ar) = 0.0092;
    c(N2) = 1 - c(O2) - c(Ar);
end
if(N == 2)
    % Use binary air.
    c(O2) = 0.21;
    c(N2) = 1 - c(O2);
end

% Load a file to save recalculating things up to here.
load Tsx_Data

% We will build the surface as we go.  Start a figure and put it on hold.
figure(1)
clf
hold on
ylabel('Nitrogen Mole Fraction','rotation',0)
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
zlabel('Temperature (K)')
view([-10 40])
axis([2 7 0 1 60 300])
grid on
drawnow

% Tell the user the composition planes to be shown.
xN2_planes = clist

% Show the P-rho critical locus on the plot.
plot3(sc/1e3,xN2,Tc,'kx','LineWidth',2)

% Add domes at each plane.
for i=1:1:NCS
    plot3(sdome(i,:)/1e3,xN2(i)*ones(1,length(sdome)),Tdome(i,:),'k-','LineWidth',2)
    drawnow
end

% Add dew and bubble points.
for i=1:1:NCS
    plot3(sdew(i,1)/1e3,xN2(i)*ones(1,NPSSC),Tdew(i,1),'b.','LineWidth',2)
    plot3(sbub(i,1)/1e3,xN2(i)*ones(1,NPSSC),Tbub(i,1),'b.','LineWidth',2)
    plot3(sdew(i,3)/1e3,xN2(i)*ones(1,NPSSC),Tdew(i,3),'b.','LineWidth',2)
    plot3(sbub(i,3)/1e3,xN2(i)*ones(1,NPSSC),Tbub(i,3),'b.','LineWidth',2)
    drawnow
end

% Put splines through the bubble and dew data in the composition direction.
[row col] = size(Tdew);
for j=1:1:col
    plot3(sdewSpline(1,:)/1e3,xN2Spline,TdewSpline(1,:),'b-','LineWidth',2)
    plot3(sbubSpline(1,:)/1e3,xN2Spline,TbubSpline(1,:),'b-','LineWidth',2)
    plot3(sdewSpline(3,:)/1e3,xN2Spline,TdewSpline(3,:),'b-','LineWidth',2)
    plot3(sbubSpline(3,:)/1e3,xN2Spline,TbubSpline(3,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist

% Add isobars.
for i=1:1:NCS
    % Do the points below the dome.
    for j=1:2:3
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
%     % Now do the pressures that do not have dew points (above dome).
%     for j=j+1:1:NPS
%         plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
%         plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
%         drawnow
%     end
end

% Finished with the dynamic plot.  Close it out.
%hold off