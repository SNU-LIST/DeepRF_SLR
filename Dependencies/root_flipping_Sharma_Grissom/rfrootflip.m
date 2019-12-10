function [rfout,bestroots,bestflips,initroots] = rfrootflip(b,ntrials,bc,flip,d1,TBW)

% Beta polynomial root flipping 
% Inputs:
%   b: beta coefficients
%   ntrials: # monte-carlo trials
%   bc: normalized band centers in radians (-pi to pi)
%   flip: flip angle in radians
%   d1: passband ripple
%   TBW: time-bandwidth product of pulse
% Outputs:
%   rfout: rf with lowest peak amplitude
%   bestroots: roots corresponding to rfout
%   bestflips: flipped root pattern
%   initroots: initial roots
% 
% Adapted from Miki Lustig's flipZeros code by Will Grissom and Anuj Sharma
% Vanderbilt University, 2015

N = length(b); % number of time points in pulse

b = b./max(abs(ft(b))); % normalized beta coefficients
b = b*sin(flip/2 + atan(d1*2)/2); % scale to target flip angle
rfinit_max = max(abs(b2rf(b))); % get initial peak RF amplitude
r = roots(b); % calculate roots of beta polynomial

% sort the roots by phase angle
[~,idx] = sort(angle(r)); 
r = r(idx);

iter_accum = 0;
convergence = 10000;

minpeakrf = Inf;
for ii = 1:ntrials
    
    % determine which indices to flip
    doflip = flipZeros(N-1,r,(TBW/N*pi)*3,bc-pi/N);
   
    % flip those indices
    rt = r;
    rt(doflip == 1) = conj(1./rt(doflip==1));
    
    % get root-flipped RF
    R = poly(leja_fast(rt)); % get polynomial coefficients back from flipped roots
    R = R/max(abs(freqz(R))); % normalized by max of filter response
    bt = R*sin(flip/2 + atan(d1*2)/2); % scale to target flip
    rft = b2rf(bt); % calculate the RF pulse from root-flipped beta

    % save result if better
    if max(abs(rft)) < minpeakrf
        minpeakrf = max(abs(rft));
        bestflips = doflip;
        rfout = rft;
        bestroots = rt;
        iter_accum = 0;
    end
    
    if rem(ii,500) == 0
        fprintf('Iteration %d of %d. Peak RF: %0.2d rad. Peak init RF: %0.2d rad.\n', ...
            ii,ntrials,minpeakrf,rfinit_max);
    end
    
    iter_accum = iter_accum + 1;
    
    if iter_accum > convergence
        disp('end of execution');
        break;        
    end
    
end

initroots = r;