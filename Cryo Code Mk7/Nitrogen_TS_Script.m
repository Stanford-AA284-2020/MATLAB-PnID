% Temperature-entropy diagram for Nitrogen

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

% Make the plot.
figure(1)
clf
hold on

fprintf('Plotting T-s curves...\n')

% Get the critical and triple data.
Tc = Tcrit_i(N2);
rc = rcrit_i(N2);
Pc = Pcrit_i(N2);
sc = s_irT(N2,rc,Tc);
Tt = Ttrip_i(N2);
Pt = Ptrip_i(N2);
sft = s_irT(N2,rftrip_i(N2),Tt);
sgt = s_irT(N2,rgtrip_i(N2),Tt);

% Get the numerical critical and triple data.
Tc_num = 126.192;
rc_num = 313.300;
Pc_num = P_irT(N2,rc_num,Tc_num);
sc_num = s_irT(N2,rc_num,Tc_num);
[Pt_num rft_num rgt_num] = Saturation_iT(N2,Tt);
sft_num = s_irT(N2,rft_num,Tt);
sgt_num = s_irT(N2,rgt_num,Tt);

% Set the minimum T props.
Tmin = Tt;
Pmin = Pt_num;

% Set the limits for data curves.
Tmax = 600;
Pmax = 1000*100000;

% Make a vapor dome.
dT = (Tc_num-Tmin)/100;
i = 1;
% Preallocate storage...
Tsatline = zeros(101);
sliqline = zeros(101);
svapline = zeros(101);
for T=Tmin:dT:Tc_num-dT    % Stop short of the critical point.
    Tsatline(i) = T;
    [Psat rf rg] = Saturation_iT(N2,T);
    sliqline(i) = s_irT(N2,rf,T);
    svapline(i) = s_irT(N2,rg,T);
    i = i+1;
end
Tsatline(i) = Tc;   % Add the critical point now.
sliqline(i) = sc;
svapline(i) = sc;

% Start a set of isobaric curves.
Plist = [0.5 1 2 5 10 20 50 70 100 200 500 1000];  % Pressure in bar
% Preallocate storage...
Tpresline = zeros(length(Plist),50);
spresline = zeros(length(Plist),50);
for j=1:1:6
    P = Plist(j)*100000

    % Find the saturation state at this pressure.
    [Tsat rf rg] = Saturation_iP(N2,P); 
    sf = s_irT(N2,rf,Tsat);
    sg = s_irT(N2,rg,Tsat);
    
    % Do the compressed liquid side.
    dT = (Tsat - Tmin)/50;
    i = 1;
    for T=Tmin:dT:Tsat-dT  % Stop short of saturation.
        Tpresline(j,i) = T;
        r = rl_iTP(N2,T,P);
        spresline(j,i) = s_irT(N2,r,T);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    spresline(j,i) = sf;
    i = i+1;

    % Now go across the dome.
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Tpresline(j,i) = Tsat;
        spresline(j,i) = sf + q*(sg-sf);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    spresline(j,i) = sg;
    i = i+1;

    % Do the vapor side.
    dT = (Tmax - Tsat)/150;
    for T=Tsat+dT:dT:Tmax  % Start just above saturation.
        Tpresline(j,i) = T;
        r = rv_iTP(N2,T,P);
        spresline(j,i) = s_irT(N2,r,T);
        i = i+1;
    end
end

% Add isobars above the critical pressure.
dT = (Tmax - Tmin)/250;
for j=j+1:1:length(Plist)
    P = Plist(j)*100000;       % In Pascals
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tpresline(j,i) = T;
        r = rl_iTP(N2,T,P);
        spresline(j,i) = s_irT(N2,r,T);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add some isenthalpic curves.
hsteps = 15;
Psteps = 50;
r_min = rl_iTP(N2,80,0.5E05);
hlow = h_irT(N2,r_min,80);
r_high = rv_iTP(N2,600,1E05);
hhigh = h_irT(N2,r_high,600);
dh = (hhigh - hlow)/hsteps;
Phigh = Plist(length(Plist))*1e5;
Plow  = Plist(1)*1e5;
Plist_plot = logspace(log10(Phigh), log10(Plow), Psteps);
j = 1;
% Preallocate storage...
hlist       = zeros(hsteps+1);
Tenthline   = zeros(hsteps+1,Psteps);
senthline   = zeros(hsteps+1,Psteps);
for h=hhigh:-dh:hlow
    h
    hlist(j) = h/1e6;
    for i=1:1:Psteps
        P = Plist_plot(i);
        % Watch out for enthalpy below Tmin
        rmin = rl_iTP(N2,Tmin,P);
        hmin = h_irT(N2,rmin,Tmin);
        smin = s_irT(N2,rmin,Tmin);
        if(h > hmin)
            [r T rf rg] = rT_ihP(N2,h,P);
            if(rf == 0) % Not under the dome.
                Tenthline(j,i) = T;
                senthline(j,i) = s_irT(N2,r,T);
            else        % Under the dome.
                Tenthline(j,i) = T;
                q = (1/r-1/rf)/(1/rg-1/rf);
                sf = s_irT(N2,rf,T);
                sg = s_irT(N2,rg,T);
                senthline(j,i) = sf + q*(sg-sf);
            end
        else
            Tenthline(j,i) = Tmin;
            senthline(j,i) = smin;
        end
    end
    j = j+1;
end
N_lines = hsteps+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the plot.
figure(1)
clf
hold on
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,Tpresline(i,:),'Color',[.3 .3 1])
end
for i=1:1:N_lines
    plot(senthline(i,:)/1000,Tenthline(i,:),'Color',[0 .6 0])
end
plot(sliqline/1000,Tsatline,'k')
plot(svapline/1000,Tsatline,'k')
plot([sft_num sgt_num]/1000, [Tt Tt],'k')
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Temperature (K)')
for i=1:1:length(Plist)-1
    text(spresline(i,251)/1000,Tpresline(i,251)+6,...
        num2str(Plist(i)),'Color',[.1 .1 1])
end
i = length(Plist);
text(spresline(i,251)/1000-6,Tpresline(i,251)+6,...
    ['P = ',num2str(Plist(i)),' bar'],'Color',[.1 .1 1])
for i=1:1:N_lines-1
    text(senthline(i,Psteps)/1000,Tenthline(i,Psteps)-2,...
        num2str(hlist(i)),'Color',[0 .5 0])
end
i = length(hlist);
text(senthline(i,Psteps)/1000,Tenthline(i,Psteps)-2,...
    ['h = ',num2str(hlist(i)),' MJ/kg'],'Color',[0 .5 0])
hold off
plotfixer
