function rf90 = matched_exc_design(rf180)

% function to compute matched excitation pulse for spin-echo imaging
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


N = length(rf180); % # time points
[~,b90di] = rf2ab_sas(rf180,(-N/2:1/2:N/2-1/2)',0); % get beta of 180 pulse
b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
b90d = -conj(b90d);
bx=fftshift( fft(ifftshift(b90d))/length(b90d) );
[~,ax] = b2amp(bx);
rf90 = -1i*conj(islr(ax,bx));
rf90=real(rf90);
