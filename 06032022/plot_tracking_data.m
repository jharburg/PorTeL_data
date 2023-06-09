%% analyze 06/03/22 ISS pass sim

%% load data
%simulation track (.1 sec step, from ref mount pos)
% track_TLE = readmatrix("TLE_track_20220603_131404.38.csv");

%simulation track (.1 sec step, from prev mount pos)
% track_TLE = readmatrix("TLE_track_20220603_153840.03.csv");

%simulation track (.2 sec step, from prev mount pos)
% track_TLE = readmatrix("TLE_track_20220603_155502.86.csv");

%simulation track (.2 sec step, from prev mount pos, removed .15sec offset)
% track_TLE = readmatrix("TLE_track_20220603_162923.30.csv");

%simulation track (timed step, from prev mount pos, removed .15sec offset)
% track_TLE = readmatrix("TLE_track_20220613_162723.23.csv");

%simulation track (timed step, from prev mount pos, removed .15sec offset, changed gains to 0.2)
% track_TLE = readmatrix("TLE_track_20220614_102548.70.csv");

%simulation track (timed step (same time as sat vec calc), from prev mount pos, added .15sec offset, changed gains to 0.2)
% track_TLE = readmatrix("TLE_track_20220628_151833.40.csv");

%simulation track (timed step (same time as sat vec calc), from prev mount pos, removed .15sec offset, changed gains to 0.2)
% track_TLE = readmatrix("TLE_track_20220628_145916.03.csv");

%simulation track (timed step (same time as sat vec calc), from prev mount pos, removed .15sec offset, changed gains to 0.2, added I term)
track_TLE = readmatrix("TLE_track_20220628_155414.83.csv");

%% plot TLE
t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

%  az/al vs time
figure(1)
plot(t,track_TLE(:,2:5))
xlabel('Time (sec)')
ylabel('Gimbal Angles (deg)')
legend('measured az','measured el','reference az','reference el')

%  az/al error vs time
figure(2)
plot(t,(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')

%  az/al rate vs time
figure(3)
plot(t,track_TLE(:,6:7)/3600)
hold on
plot(t,track_TLE(:,8:9))
hold off
% plot(t,track_TLE(:,6:7)/3600-track_TLE(:,8:9))
xlabel('Time (sec)')
ylabel('Gimbal Rates (deg/sec)')
legend('commanded az rate','commanded el rate','reference az rate','reference el rate')

% show time gaps in control step
figure(4)
plot(t(1:end-1),diff(t))
xlabel('Time (sec)')
ylabel('control step dt (sec)')

%% test contol methods
%simulation track (timed step (same time as sat vec calc), from prev mount pos, removed .15sec offset, changed gains to 0.2)
track_TLE = readmatrix("TLE_track_20220628_145916.03.csv");

t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = ref_ang(1,:);
dt = diff(t);
dv_p = 0;
dv_i = 0;

% ref_ang_off = ref_ang; %test target offset
% ns = 5; %# steps offset (~.2 sec steps)
% ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
% dt_off = diff(t);
% dt_off(1:end-ns,:) = dt_off(ns+1:end,:);

for i=2:length(ref_rate)

    %reference rate
%     v = 0;
    v = ref_rate(i-1,:);
%     v = (ref_rate(i-1,:)+ref_rate(i,:))/2;

    %correction (proportional term)
%     ref_ang = ref_ang_off;  %test target offset
%     dt = dt_off;  %test target offset
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
%     if max(abs(dv_p))<=60/3600
%         dv_i = dv_i + dv_p;
%     end

    %update step
%     k=[0.6,1];
    k=[0.2,0.2];
%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*(v+1./k.*(dv_p+1*dv_i));

end

% %  az/al vs time
% figure()
% plot(t,[sim_ang ref_ang])
% xlabel('Time (sec)')
% ylabel('Gimbal Angles (deg)')
% legend('sim az','sim el','reference az','reference el')

figure()
plot(t,(sim_ang-ref_ang)*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
