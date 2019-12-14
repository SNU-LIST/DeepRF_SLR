function [root_num,idxPass,r] = num_root(bc,TBW,r,N)

% function to count the number of eligible roots
%
% Subsantial amount of code adopted from here:
%   1) http://www.vuiis.vanderbilt.edu/~grissowa/
%   2) http://rsl.stanford.edu/research/software.html
%   3) https://github.com/mriphysics/AM_multiband/
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

N = N-1; % N is always 512
wp =(TBW/N*pi)*3;
bp = bc-pi/N;

% find the indices of the bands
w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
end

% eligible roots for flipping are on the top half of complex plane
index = idxPass(1:find(idxPass == (N+1)/2));
root_num = size(index,2);

% visualize complex plane
% figure; scatter(real(r),imag(r)); hold on;
% scatter(real(r(index)),imag(r(index)));
% legend('roots','eligible roots');
% pause(0.1); pause;

end