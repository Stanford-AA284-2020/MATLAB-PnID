% Temperature-entropy diagram for Oxygen


% Provide access to support files via the Matlab path.
addpath 'Fundamental Relation Files' 
addpath 'Fundamental Relation Data'
addpath 'Setup Files' 
addpath 'Property Files' 

% Clean up and get ready to go.
clear all; close all; clc;
format compact
fprintf('\n**************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Make the plot.
figure(1)
clf
hold on

fprintf('Plotting P-v curves...\n')

% Get the critical and triple data.
Tc = Tcrit_i(O2);
rc = rcrit_i(O2);
Pc = Pcrit_i(O2);
vc = 1/rc;
Tt = Ttrip_i(O2);
Pt = Ptrip_i(O2);
rftrip = rftrip_i(O2);
rgtrip = rgtrip_i(O2);

% Get the numerical critical and triple data.
Tc_num = 154.5994;
rc_num = 426.9340; 
vc_num = 1/rc_num;
Pc_num = P_irT(O2,rc_num,Tc_num);
sc_num = s_irT(O2,rc_num,Tc_num);
[Pt_num rft_num rgt_num] = Saturation_iT(O2,Tt);
vftrip = 1/rft_num; 
vgtrip = 1/rgt_num;

% Set the minimum T props.
Tmin = Tt;
Pmin = 0.01*1E5;

% Set the limits for data curves.
Tmax = 1000;
Pmax = 1000*100000;

% Make a vapor dome.
dP = (Pc_num-Pmin)/500;
i = 1;
% Preallocate storage...
Psatline = zeros(501);
vliqline = zeros(501);
vvapline = zeros(501);
for P =Pmin:dP:Pc_num-dP    % Stop short of the critical point.
    i
    Psatline(i) = P;
    [Tsat rf rg] = Saturation_iP(O2,P);
    vliqline(i) = 1/rf;
    vvapline(i) = 1/rg;
    i = i+1;
end
Psatline(i) = Pc;   % Add the critical point now.
vliqline(i) = 1/rc_num;
vvapline(i) = 1/rc_num; % Fine until here

% Start a set of isotherms
Tlist = [90 100 110 120 130 140 180 200 250 300 500 800 1000];
% Preallocate storage...
Ppresline = zeros(length(Tlist),300);
vpresline = zeros(length(Tlist),300);

for j=1:1:6
    
    T = Tlist(j);
    % Find the saturation state at this pressure.
    [Psat rf rg] = Saturation_iT(O2,T); 
   
   % Do the compressed liquid side.
    dP = (Pmax - Psat)/100;
    i = 1;
    for P=Pmax:-dP:Psat+dP  % Stop short of saturation.
        Ppresline(j,i) = P;
        r = rl_iTP(O2,T,P);
        vpresline(j,i) = 1/r;
        i = i+1;
    end
    Ppresline(j,i) = Psat;   % Add the saturation point now.
    vpresline(j,i) = 1/rf;
    i = i+1;

    % Now go across the dome.
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Ppresline(j,i) = Psat;
        vpresline(j,i) = (q/rg) + ((1-q)/rf);
        i = i+1;
    end
    Ppresline(j,i) = Psat;   % Add the saturation point now.
    vpresline(j,i) = 1/rg;
    i = i+1;

    % Do the vapor side.
    dP = (Psat - Pmin)/150;
    for P=Psat-dP:-dP:Pmin  % Start just below saturation.
        Ppresline(j,i) = P;
        r = rv_iTP(O2,T,P);
        vpresline(j,i) = 1/r;
        i = i+1;
    end
end

% Add isotherms above the critical temperature.
dP = (Pmax - Pmin)/301;
for j=j+1:1:length(Tlist)
    T = Tlist(j);
    i = 1;
    for P=Pmax:-dP:Pmin  % Stop short of saturation.
        Ppresline(j,i) = P;
        r = rl_iTP(O2,T,P);
        vpresline(j,i) = 1/r;
        i = i+1;
    end
end
Ppresline = Ppresline(:, 1:300);
vpresline = vpresline(:, 1:300);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the plot.
figure(1)
clf
semilogx(vliqline,Psatline/1E5,'k');
hold on;
semilogx(vvapline,Psatline/1E5,'k')
hold on

for i=1:1:length(Tlist)
    semilogx(vpresline(i,:),Ppresline(i,:)/1E5,'r-')
end
axis([5E-4 0.5 0.05 200]);
xlabel('Specific Volume [m^3/kg]'); ylabel('Pressure [bar]');
title('Process Plot (PV)');
plotfixer;

save('Oxygen_PV.mat', 'vliqline', 'vvapline', 'Psatline', 'vpresline', 'Ppresline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Add some isenthalpic curves.
% hsteps = 50;
% r_min = rl_iTP(O2,110,1000E05);
% hlow = h_irT(O2,r_min,110);
% r_high = rv_iTP(O2,600,1E05);
% hhigh = h_irT(O2,r_high,600);
% dh = (hhigh - hlow)/hsteps;
% Plist = [0.1 0.5 1 2 5 10 20 50 70 100 200 500 1000]*1E5;  
% j = 1;
% % Preallocate storage...
% hlist       = zeros(hsteps+1);
% Penthline   = zeros(hsteps+1,length(Plist));
% venthline   = zeros(hsteps+1,length(Plist));
% for h=hhigh:-dh:hlow
%     h
%     hlist(j) = h;
%     for i=1:1:length(Plist)
%         P = Plist(i);
%         % Watch out for enthalpy below Tmin
%         rmin = rl_iTP(O2,Tmin,P);
%         hmin = h_irT(O2,rmin,Tmin);
%         smin = s_irT(O2,rmin,Tmin);
%         if(h > hmin)
%             [r P rf rg] = rT_ihP(O2,h,P);
%             if(rf == 0) % Not under the dome.
%                 Penthline(j,i) = P;
%                 venthline(j,i) = 1/r;
%             else        % Under the dome.
%                 Penthline(j,i) = P;
%                 
%                 venthline(j,i) = r;
%             end
%         else
%             Penthline(j,i) = Pmin;
%             venthline(j,i) = 1/rmin;
%         end
%     end
%     j = j+1;
% end
% N_lines = hsteps+1;
% 
% figure(1);
% hold on;
% for i=1:1:N_lines
%    semilogx(venthline(i,:),Penthline(i,:),'Color',[0 .6 0]);hold on;
% end
% 
% % plot([sft_num sgt_num]/1000, [Tt Tt],'k')
% % xlabel('Specific Entropy (kJ/kg-K)')
% % ylabel('Temperature (K)')
% % for i=1:1:length(Plist)-1
% %     text(vpresline(i,251)/1000,Ppresline(i,251)+6,...
% %         num2str(Plist(i)),'Color',[.1 .1 1])
% % end
% % i = length(Plist);
% % text(vpresline(i,251)/1000-6,Ppresline(i,251)+6,...
% %     ['P = ',num2str(Plist(i)),' bar'],'Color',[.1 .1 1])
% % for i=1:1:N_lines-1
% %     text(senthline(i,Psteps)/1000,Tenthline(i,Psteps)-2,...
% %         num2str(hlist(i)),'Color',[0 .5 0])
% % end
% % i = length(hlist);
% % text(senthline(i,Psteps)/1000,Tenthline(i,Psteps)-2,...
% %     ['h = ',num2str(hlist(i)),' MJ/kg'],'Color',[0 .5 0])
% % hold off
% % plotfixer
