function [a,b] = blochsim(rf,gamgdt,xx)

[Ns,Nd] = size(xx); % Ns: # spatial locs. Nd: # of dimensions
Nt = size(gamgdt,1); % Nt: # time points
a = ones(Ns,1);
b = zeros(Ns,1);

for ii = 1:Nt
    
  % apply rf
  C = cos(abs(rf(ii))/2);
  S = 1i*exp(1i*angle(rf(ii)))*sin(abs(rf(ii))/2);
  at = a*C - b*conj(S);
  bt = a*S + b*C;

  a = at;b = bt;
  
  % apply phase accrual
  z = exp(-1i*xx*gamgdt(ii,:)');  
  b = b.*z;
  
end

z = exp(1i/2*xx*sum(gamgdt,1)');
a = a.*z;
b = b.*z;

