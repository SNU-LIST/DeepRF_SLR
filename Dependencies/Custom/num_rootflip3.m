function [root_num,idxPass,r] = num_rootflip3(bc,TBW,r,N)

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