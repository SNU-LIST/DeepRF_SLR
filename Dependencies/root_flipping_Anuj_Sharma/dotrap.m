function [trap,ramppts]=dotrap(area,gmax,dgdt,dt);
% trapezoid for ramp-sampled RF
% dotrap(area,gmax,dgdt,dt);
%   area = pulse area in (g sec)/cm
%   gmax = max gradient in g/cm
%   dgdt = max slew in g/cm/sec
%   dt   = sample time in sec

if abs(area) > 0

  ramppts=ceil(gmax/dgdt/dt);
  triareamax = ramppts*dt*gmax;

  if triareamax >= abs(area)
    % triangle pulse
    newgmax = sqrt(abs(area)*dgdt);
    ramppts=ceil(newgmax/dgdt/dt);
    pulse = [(0:ramppts)./ramppts (ramppts:-1:0)./ramppts];
  else
    % trap pulse
    nflat = ceil((area-triareamax)/gmax/dt);
    pulse = [(0:ramppts)./ramppts ones([1 nflat]) (ramppts:-1:0)./ramppts];
  end
  
  trap = pulse*(area/(sum(pulse)*dt));
  
else
  trap = 0;
end

