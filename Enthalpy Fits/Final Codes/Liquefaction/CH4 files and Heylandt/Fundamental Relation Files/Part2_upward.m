%% Part 2 upward
close all
clear all

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

N = 3;
Setup_Air_Props;

% load data: output state = bottom product + vapor out of reboiler and into tray
% load data: input state = liquid going into reboiler and output vapor from tray

% load part1.mat  % load inlet(vapor leaving) and outlet (liquid coming down)
load part1reboiler.mat % load vapor_in (vapor entering) and liquid_out (liquid leaving below)

below = state_out;
tray = state_in;   % vapor entering this tray

qout = 0:0.1:1;
q_tray = zeros(length(qout),1);
above(length(qout))=struct('T',[],'P',[],'x',[],'y',[],'rf',[],'rg',[],'V',[],'L',[]);

for i = 1:length(qout)
    % using reboiler quality to calculate flow rate entering from below
    if qout(i) == 1
        % when quality is 1, no liquid coming out
        below.L = 0;
        below.V = 1;
        tray(i).V = below.V; % vapor pass through, no heat transfer
    else
        % when quality is not 1, set liquid exiting the bottom to 1 and find
        % others
        below.L = 1;
    end
    
    % set starting point for tray quality
    if i == 1
        q_tray_guess = 0.5;
    elseif i == length(qout)
        q_tray_guess = 1;
    else
        q_tray_guess = q_tray(i-1);
    end
    
    [tray(i), above(i),q_tray(i)] = upwd_tray(above,tray(i),below(i),q_tray_guess);
    
end 

save('part2reboiler.mat')


