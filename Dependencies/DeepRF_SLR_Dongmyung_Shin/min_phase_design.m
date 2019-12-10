function [bp,d1,tb,r,N] = min_phase_design(n,tblin,d1,d2,bandsep,nb)

% function to design initial minimum phase pulse for root-flipping
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

% directly design a MB min-phase filter, whose roots we can flip
d1 = d1/4; % target beta passband ripple, considering full 90-180 pair
d2 = (d2/sqrt(2))^0.25; % target beta stopband ripple, considering full 90-180 pair

dinfmin = 1/2*dinf(2*d1,d2^2/2); % d-infinity for a min-phase pulse with beta ripples (d1,d2)
dinflin = dinf(d1,d2);      % d-infinity for a linear phase pulse with the same ripples
tb = tblin/dinflin*dinfmin; % scale TBW product so as to get the same transition
% width as linear phase pulse with same ripples,
% after scaling back to desired slice thickness. This
% makes comparison to other MB excitations more
% meaningful, since all will have same slice characteristics.
w = dinfmin/tb; % transition width

N = 2*(n-1); % length of linear-phase filter we will factor to get min-phase filter

% use cvx as described in the MRM paper
osfact = 10; % oversampling factor

nn = (0:N/2*osfact)'/(N/2*osfact);  % 0 to 1 - profile indices
d = zeros(N/2*osfact+1,1);          % passband mask
s = zeros(N/2*osfact+1,1);          % stopband mask
wts = zeros(N/2*osfact+1,1);        % ripple taper weights

if rem(nb,2) % if Nb odd
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
    d = nn <= f(2)/(n/2); % target pattern
    wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
else
    f = 0;
end

% add non-DC bands to the profiles
for ii = 1:floor(nb/2)
    cent = (ii - (rem(nb,2) == 0)/2)*(bandsep)*dinfmin/dinflin;
    f = [f (cent-(1+w)*(tb/2)) (cent-(1-w)*(tb/2)) (cent+(1-w)*(tb/2)) (cent+(1+w)*(tb/2))];
    d = d | (nn >= f(end-2)/(n/2) & nn <= f(end-1)/(n/2));
    s = s | (nn >= f(end-4)/(n/2) & nn <= f(end-3)/(n/2));
    nnc = nn - (f(end-1)+f(end-2))/2/(n/2); % indices centered with passband for weight calcs
    wts = max(wts,1./abs(nnc).^2); % quadratically-decaying ripple weights
end
% append the last stopband
s = s | (nn >= f(end)/(n/2));
wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

% build system matrix for cvx design
A = 2*cos(2*pi*(0:N/2*osfact)'*(-N/2:0)/(N*osfact));A(:,end) = 1/2*A(:,end);

% mask everything to get rid of transition bands
Ad = A(s | d,:);
dd = double(d(s | d));
ss = wts(s | d).*double(s(s | d));

% use cvx to do the constrained optimization
cvx_begin
variable delta(1)
variable x(N/2+1)
minimize( delta )
subject to
-delta*dd <= Ad*x - dd <= delta*dd + delta*d2^2/(2*d1)*ss
cvx_end

% stack two halves together to get full linear-phase filter
x = [x;x(end-1:-1:1)]';


% factor the linear phase filter to get a min-phase filter h
h = real(fmp(x));
h = h(end:-1:1);

% root flip h!
bp = 2*((1:nb)-1/2-nb/2)*bandsep/n*pi*dinfmin/dinflin; % centers of passbands

b = h;

N = length(b); % number of time points in pulse

b = b./max(abs(ft(b))); % normalized beta coefficients
b = b*sin(pi/2 + atan(d1*2)/2); % scale to target flip angle
r = roots(b); % calculate roots of beta polynomial

% sort the roots by phase angle
[~,idx] = sort(angle(r));
r = r(idx);

end

