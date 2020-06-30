% origianlly from Minju Jo
% modified by Dongmyung Shin
% Seoul National University
% Laboratory for Imaging Science and Technology (LIST), 2019


function [Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,RF,Gz,time_step,slice_thick,spatial_point)

% units
% x, y, z : (spatial_point*1) for each
% T1, T2 : [msec]
% RF (2*N) : RF(1,:) is [Gauss], RF(2,:) row is angle [degree]
% Gz : [mT/m] (1 * N matrix) = [0.1 G/cm]
% time_step : pulse duration / sample points [msec]
% slice_thick : half of simulating thickness [m]
% spatial_point : number of points that are going to be simulated [#]
% 
% returns
% Mx, My, Mz : (off resonance * time series)

% use units in standard unit (sec, Hz, etc.)
T1=T1/1000; % [sec]
T2=T2/1000; % [sec]
length_RF = size(RF,2); % 2(mag, phase)
t_int = time_step*10^(-3); % [sec]
delta_omega = 2*pi*42.57747892*10^6 * (-1:2/(spatial_point-1):1).'*slice_thick*0.001*Gz;

Mx = repmat(x,[1 length_RF]);
My = repmat(y,[1 length_RF]); 
Mz = repmat(z,[1 length_RF]); 

RF_amp = repmat(RF(1,:)*2*pi*4257.747892,[spatial_point 1]);
RF_phase = repmat(RF(2,:),[spatial_point 1])*pi/180;

alpha = t_int*sqrt(RF_amp.^2+delta_omega.^2);
zeta = atan2(RF_amp, delta_omega);
theta = RF_phase;

ca = cos(alpha); sa = sin(alpha);
cz = cos(zeta); sz = sin(zeta);
ct = cos(theta); st = sin(theta);
E1=exp(-t_int/T1); E2=exp(-t_int/T2);

Mx_x_part = ct.*(E2.*ct.*sz.^2 + cz.*(E2.*sa.*st + E2.*ca.*ct.*cz)) + st.*(E2.*ca.*st - E2.*ct.*cz.*sa);
Mx_y_part = st.*(E2.*ct.*sz.^2 + cz.*(E2.*sa.*st + E2.*ca.*ct.*cz)) - ct.*(E2.*ca.*st - E2.*ct.*cz.*sa);
Mx_z_part = E2.*ct.*cz.*sz - sz.*(E2.*sa.*st + E2.*ca.*ct.*cz);
My_x_part = - ct.*(- E2.*st.*sz.^2 + cz.*(E2.*ct.*sa - E2.*ca.*cz.*st)) - st.*(E2.*ca.*ct + E2.*cz.*sa.*st);
My_y_part = ct.*(E2.*ca.*ct + E2.*cz.*sa.*st) - st.*(- E2.*st.*sz.^2 + cz.*(E2.*ct.*sa - E2.*ca.*cz.*st));
My_z_part = sz.*(E2.*ct.*sa - E2.*ca.*cz.*st) + E2.*cz.*st.*sz;
Mz_x_part = ct.*(E1.*cz.*sz - E1.*ca.*cz.*sz) + E1.*sa.*st.*sz;
Mz_y_part = st.*(E1.*cz.*sz - E1.*ca.*cz.*sz) - E1.*ct.*sa.*sz;
Mz_z_part = E1.*cz.^2 + E1.*ca.*sz.^2;

for q=1:length_RF-1
    Mx(:,q+1) = Mx(:,q).*Mx_x_part(:,q) + My(:,q).*Mx_y_part(:,q) + Mz(:,q).*Mx_z_part(:,q);
    My(:,q+1) = Mx(:,q).*My_x_part(:,q) + My(:,q).*My_y_part(:,q) + Mz(:,q).*My_z_part(:,q);
    Mz(:,q+1) = Mx(:,q).*Mz_x_part(:,q) + My(:,q).*Mz_y_part(:,q) + Mz(:,q).*Mz_z_part(:,q) + (1-E1);
end

end