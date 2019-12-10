% This script is to design multiband refocusing pulse
% using Monte-Carlo algorithm
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


%%
% design parameters
gamma = 2*pi*4258;      % 1H gyromagnetic ratio (radians/second/gauss)
n = 512;                % # of time points
nb = 3;                 % # of bands
tb = 6;                 % Time-bandwidth product
slthick = 5;            % mm, slice thickness for simulation
bandsep = 6*tb;         % Slice gap in unit of slice thickness
d1e = 0.01;             % combined Mxy ripple, passband
d2e = 0.01;             % combined Mxy ripple, stopband
maxb1 = 0.20;           % Gauss, system b1 limit
ntrials = 1000;   % flipping trials

% initial design of minimum phase SLR pulse
[bp,d1,tbrf,rr,N] = min_phase_design(n,tb,d1e,d2e,bandsep,nb);

% how many roots can be flipped?
[root_num,idxPass,r] = num_root(bp,tbrf,rr,N);
disp(['number of eligible roots = ',num2str(root_num)]);

% tracking varialbes & arrays
maxrf = 100;
num_iteration = 0;
iter_list = []; % number of flipping trials
time_list = []; % execution time
amp_list = []; % peak RF
rf_list = []; % RF
pattern_list = []; % flipping pattern


%%
% start Monte-Carlo algorithm
tic

for i=1:ntrials
    % distribute probabilities across the pasband roots
    p = 1/length(idxPass):1/length(idxPass):1;
    flip = (rand(1,length(idxPass))<2*p);
    % acquire binary flipping pattern
    pattern = flip(1:root_num);
    % calculate RF pulse
    [rf180m,rf180a,~] = update_pulse(r,pattern,idxPass,N,d1);
    num_iteration = num_iteration + 1;

    % save result
    if maxrf > max(abs(rf180m))
        rf180 = rf180m.*exp(1i*rf180a);
        maxrf = max(abs(rf180));
        iter_list = [iter_list, num_iteration];
        time_list = [time_list, toc];
        amp_list = [amp_list, maxrf];
        rf_list = [rf_list; rf180];
        pattern_list = [pattern_list; pattern];
    end

    % print result
    if mod(i,500) == 0
        disp(['iteration ',num2str(i),': best max |RF| = ',num2str(maxrf)]);                
    end
end

time_list = [time_list, toc];


%%
% plot design results
rf180 = squeeze(rf_list(end,:));
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

save MC_design_result;
disp('end of execution');

