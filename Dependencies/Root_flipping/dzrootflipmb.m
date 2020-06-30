function [rf180,tb] = dzrootflipmb(n,tblin,d1,d2,bandsep,nb,type,filtAlg)

% Function to design a root-flipped multiband refocusing pulse
%
% Inputs:
%   n           # time points in the pulse waveform
%   tblin       equivalent linear phase time-bandwidth product for each band
%   (d1,d2)     maximum pass- (d1) and stop- (d2) band ripples
%               (not strictly guaranteed to be met)
%   bandsep     normalized slice separation (integer)
%   nb          number of bands
%   type        'matched' or '2xrefd'; switch to use ripple relationships for
%               matched 90-180, or twice-refocused 180
% Outputs:
%   rf180       final RF pulse in units of radians
%
% Copyright 2015, Will Grissom and Anuj Sharma, Vanderbilt University

% Check arguments

% default RF parameters
if ~exist('n','var')
    disp 'No arguments supplied; using defaults'
    n = 256;              % number of time points
    tblin = 4;            % time-bandwidth product (linear-phase)
    d1 = 0.01;            % passband Mxy ripple
    d2 = 0.01;            % stopband Mxy ripple
    bandsep = 4*tblin;    % slice gap
    nb = 6;               % number of bands
    type = '2xrefd';      % 'matched' or '2xrefd'
end
if ~exist('filtAlg','var')
    filtAlg = 'firls'; % 'firls' or 'cvx'
end


% directly design a MB min-phase filter, whose roots we can flip
if strcmp(type,'matched')
    d1 = d1/4; % target beta passband ripple, considering full 90-180 pair
    d2 = (d2/sqrt(2))^0.25; % target beta stopband ripple, considering full 90-180 pair
else
    d1 = d1/8; % target beta passband ripple, considering twice-refocused
    d2 = d2.^(1/4); % target beta stopband ripple, considering twice-refocused
end

dinfmin = 1/2*dinf(2*d1,d2^2/2); % d-infinity for a min-phase pulse with beta ripples (d1,d2)
dinflin = dinf(d1,d2);      % d-infinity for a linear phase pulse with the same ripples
tb = tblin/dinflin*dinfmin; % scale TBW product so as to get the same transition
                            % width as linear phase pulse with same ripples,
                            % after scaling back to desired slice thickness. This
                            % makes comparison to other MB excitations more
                            % meaningful, since all will have same slice characteristics.
w = dinfmin/tb; % transition width

N = 2*(n-1); % length of linear-phase filter we will factor to get min-phase filter

if strcmp(filtAlg,'firls')
    
    % use MATLAB's firls (much faster, but tends to yield a bit higher power)
    
    % build vectors of band edges, amplitudes and weights
    centers = bandsep*dinfmin/dinflin*((0:nb-1) - (nb-1)/2); % integer center of each band
    fNoShift = [-(1+w)*(tb/2) -(1-w)*(tb/2) (1-w)*(tb/2) (1+w)*(tb/2)];
    mNoShift = [0 1 1 0];
    f = -n/2;m = 0; % no passband at the most negative frequency
    % f: normalized band edges; m: filter response at each band edge
    for kk = 1:length(centers)
        f = [f centers(kk)+fNoShift];
        m = [m mNoShift];
    end
    f = [f n/2]/(n/2); % append nyquist and normalize to (0,1)
    m = [m 0]; % append final stopband
    wts = kron(ones(1,nb),[1 d1/d2]); % band error weights
    wts = [d1/d2 wts]; % leftmost stopband
    f = f(length(f)/2:end);f(1) = 0;
    m = m(length(m)/2:end);if ~rem(nb,2);m(1) = 0;end
    wts = wts(ceil(length(wts)/2):end);
    
    % use MATLAB's firls to solve for the multiband filter
    x = firls(N,f,m,wts);
    
elseif strcmp(filtAlg,'cvx')
    
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
    
else
    error 'Unknown filter design method'
end

% factor the linear phase filter to get a min-phase filter h
h = real(fmp(x));
h = h(end:-1:1);

% root flip h!
ntrials = ceil(tb)*nb*100; % number of monte carlo trials
% ntrials = 1000000;
bp = 2*((1:nb)-1/2-nb/2)*bandsep/n*pi*dinfmin/dinflin; % centers of passbands
rf180 = rfrootflip(h,ntrials,bp,pi,d1,tb); % rf180 is the root-flipped RF
