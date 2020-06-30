function [rfmag, rfangle, rt] = calculate_flipped_refo2(rt,pattern,idxPass,N,d1)

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

