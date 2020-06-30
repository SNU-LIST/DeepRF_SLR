function rf90 = dzmatchedex(rf180)

% function rf90 = dzmatchedex(rf180)
% 
% Function to compute a matched-phase pi/2 excitation pulse from
% a pi pulse

% input:    180 RF pulse in radians
% output:   matched-phase 90 RF pulse in radians

N = length(rf180); % # time points
[~,b90di] = abr(rf180,-N/2:1/2:N/2-1/2); % get beta of 180 pulse
b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
rf90 = -1i*conj(b2rf(1/length(b90d)*ft(b90d))); % solve for 90-deg pulse
