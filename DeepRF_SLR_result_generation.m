% This script to produce the design result of DeepRF_SLR
% Run this script after running 'DeepRF_SLR.py'
%
% Subsantial amount of code adopted from here:
%   1) http://www.vuiis.vanderbilt.edu/~grissowa/
%   2) http://rsl.stanford.edu/research/software.html
%   3) https://github.com/mrphysics/AM_multiband/
%
% Corresponding paper:
%   1) A Sharma, M Lustig, and W A Grissom. 
%      Root-flipped multiband refocusing pulses. 
%      Magn Reson Med 2016.
%   2) Abo Seada, S., Price, A. N., Hajnal, J. V., & Malik, S. J. 
%      Optimized amplitude modulated multiband RF pulse design. 
%      Magn Reson Med 2017.
%
% Modified at 2018 by Dongmyung Shin, Seoul National University
% shinsae11@gmail.com, http://list.snu.ac.kr

load DeepRF_SLR_refo_design.mat

% simulation parameters
gamma = 2*pi*4258;      % 1H gyromagnetic ratio (radians/second/gauss)
n = 512;                % # of time points
maxb1 = 0.20;           % Gauss, system b1 limit
slthick = 5;            % mm, slice thickness for simulation

rf180 = squeeze(rf_list(end,1,:))';
rf90 = matched_exc_design(rf180);

dt180 = max(abs(rf180))/(gamma*maxb1);
dt90 = max(abs(rf90))/(gamma*maxb1);

T180 = dt180*n;
T90 = dt90*n*2;
fprintf('Duration of 180: %0.2f ms. Duration of 90: %0.2f ms.\n',T180*1000,T90*1000);

bw180 = tbrf/T180;
bw90 = 2*tbrf/T90;

g180 = bw180/(gamma/2/pi*slthick/10); % gauss/cm
g90 = bw90/(gamma/2/pi*slthick/10);
disp(['Excitation G = ',num2str(g90*10),' mT/m, Refocusing G =',num2str(g180*10),' mT/m']);

rf180 = (abs(rf180)./(dt180*gamma)).*exp(1i*angle(rf180));
rf90 = (abs(rf90)./(dt90*gamma)).*exp(1i*angle(rf90));

[xx,yy,zz] = simul_pulse(rf90,T90,g90,rf180,T180,g180,1001,1,0.1,3000,1);

save DeepRF_SLR_design_result.mat
