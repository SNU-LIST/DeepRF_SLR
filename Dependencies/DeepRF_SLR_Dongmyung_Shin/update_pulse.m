function [rfmag, rfangle, rt] = update_pulse(rt,pattern,idxPass,N,d1)

% function to update pulse based on the given root pattern
% (conjugate symmetry enforced)
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

idxPass = idxPass(1:find(idxPass == N/2));
pattern = logical(pattern);
% enforce conjugate symmery (i.e., amplitude modulation only)
doflip = [idxPass(pattern) N-idxPass(pattern)-1];
doflip = unique(doflip);
% do flipping
rt(doflip) = conj(1./rt(doflip));

% get root-flipped RF
R = poly(leja_fast(rt)); % get polynomial coefficients back from flipped roots
R = R/max(abs(freqz(R))); % normalized by max of filter response
bt = R*sin(pi/2 + atan(d1*2)/2); % scale to target flip
rft = b2rf(bt); % calculate the RF pulse from root-flipped beta

% return RF
rfmag = abs(rft);
rfangle = angle(rft);

end

