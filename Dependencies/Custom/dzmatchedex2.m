function rf90 = dzmatchedex2(rf180)
% Function to compute a matched-phase pi/2 excitation pulse from
% a pi pulse

% input:    180 RF pulse in radians
% output:   matched-phase 90 RF pulse in radians

N = length(rf180); % # time points

[~,b90di] = rf2ab_sas(rf180,(-N/2:1/2:N/2-1/2)',0); % get beta of 180 pulse

b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
b90d = -conj(b90d);

bx=fftshift( fft(ifftshift(b90d))/length(b90d) );

[~,ax] = b2amp(bx);

rf90 = -1i*conj(islr(ax,bx));
rf90=real(rf90);
