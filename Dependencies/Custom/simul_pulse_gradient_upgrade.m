function [xx,yy,zz,Mx,My,Mz] = simul_pulse_gradient_upgrade(rf90,T90,g90,rf180,T180,g180,spatial_point,avg_filter_size,slice_thick,cruhser_size,result_plot)

% parameters
T1 = Inf;
T2 = Inf;

%% average filter applied
b = (1/avg_filter_size)*ones(1,avg_filter_size);
a = 1;
rf90 = filter(b,a,rf90);
rf180 = filter(b,a,rf180);

%% 90 pulse
% preparation
x = zeros(spatial_point,1);
y = zeros(spatial_point,1);
z = ones(spatial_point,1);
RF1 = abs(rf90);
RF2 = angle(rf90)/pi*180; % degree
Gz = ones(size(rf90))*g90*10; % mT/m
time_gap = T90*1e+3/size(RF1,2); % ms

% simulation
[Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,[RF1;RF2],Gz,time_gap,slice_thick,spatial_point);


%% left crusher
x = Mx(:,end);
y = My(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%
% assume z-component totally disappear after right crusher
z = zeros(spatial_point,1);
%%%%%%%%%%%%%%%%%%%%%%%%

time_gap = 0.001; % 1 micro sec.
RF1 = zeros(1,cruhser_size); % 3 ms duration
RF2 = zeros(1,cruhser_size);
Gz = ones(1,cruhser_size)*35; % 35 mT/m

% simulation
[Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,[RF1;RF2],Gz,time_gap,slice_thick,spatial_point);

%% 180 pulse
% preparation
x = Mx(:,end);
y = My(:,end);
z = Mz(:,end);
RF1 = abs(rf180);
RF2 = angle(rf180)/pi*180; % degree
Gz = ones(size(rf180))*g180*10; % mT/m
time_gap = T180*1e+3/size(RF1,2); % ms

% simulation
[Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,[RF1;RF2],Gz,time_gap,slice_thick,spatial_point);

%% right crusher
x = Mx(:,end);
y = My(:,end);
z = Mz(:,end);
time_gap = 0.001; % 1 micro sec.
RF1 = zeros(1,cruhser_size); % 3 ms duration
RF2 = zeros(1,cruhser_size);
Gz = ones(1,cruhser_size)*35; % 35 mT/m

% simulation
[Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,[RF1;RF2],Gz,time_gap,slice_thick,spatial_point);


%% rephasing after 180 pulse
% preparation
x = Mx(:,end);
y = My(:,end);
z = Mz(:,end);
RF1 = zeros(1,size(rf90,2));
RF2 = zeros(1,size(rf90,2));
%%%%%%%%%%%%%%%%%%%% Gz scaled a little bit%%%%%%%%%%%%%%%%
Gz = ones(1,size(rf90,2))*(g90)*10 / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_gap = T90*1e+3/size(rf90,2);

% simulation
[Mx,My,Mz] = Bloch_simul(x,y,z,T1,T2,[RF1;RF2],Gz,time_gap,slice_thick,spatial_point);


%% plot result
if result_plot
    close all;

    x = linspace(-slice_thick*1e+2,slice_thick*1e+2,spatial_point)';

    figure;
    plot(x,abs(Mx(:,end)+1i*My(:,end))); xlabel('Position (cm)'); ylabel('|Mxy/M0|');

    figure;
    plot(x,angle(exp(1i*(angle(Mx(:,end)+1i*My(:,end))-pi)))); xlabel('Position (cm)'); ylabel('Angle (radians)');
    
%     figure;
%     plot(T90/1024*1000:T90/1024*1000:T90*1000,imag(rf90));
%     ylabel 'Gauss'
%     xlabel 'ms'
%     title(['Excitation pulse (',num2str(T90*1000),' ms)']);
%     
%     figure;
%     plot(T180/512*1000:T180/512*1000:T180*1000,real(rf180));
%     ylabel 'Gauss'
%     xlabel 'ms'
%     title(['Refocusing pulse (',num2str(T180*1000),' ms)']);

    n = size(rf180,2);

    figure;
    plot(T90/2/n*1000:T90/2/n*1000:T90*1000,rf90);
    ylabel 'Gauss'
    xlabel 'milliseconds'
    title(['Excitation pulse (',num2str(T90*1000),' ms)']);

%     figure;
%     plot(T90/2/n*1000:T90/2/n*1000:T90*1000,angle(rf90));
%     ylabel 'Radians'
%     xlabel 'milliseconds'
%     title(['Excitation pulse (',num2str(T90*1000),' ms)']);

    figure;
    plot(T180/n*1000:T180/n*1000:T180*1000,rf180);
    ylabel 'Gauss'
    xlabel 'milliseconds'
    title(['Refocusing pulse (',num2str(T180*1000),' ms)']);

%     figure;
%     plot(T180/n*1000:T180/n*1000:T180*1000,angle(rf180));
%     ylabel 'Radians'
%     xlabel 'milliseconds'
%     title(['Refocusing pulse (',num2str(T180*1000),' ms)']);
end

%% 
xx = Mx(:,end);
yy = My(:,end);
zz = Mz(:,end);

end

