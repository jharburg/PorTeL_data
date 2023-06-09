%% analyze 09/15/22

%% load data
%gains 0.25, no integral term 

%track star
track_FSM = readmatrix("FSM_CLtrack_20220915_200306.88.csv");

%track iss (wrong TLE, empty array)
% track_TLE = readmatrix("TLE_track_20220915_204513.91.csv");

%track iss (had to dq to brightest, ~1 deg off)
% track_TLE = readmatrix("TLE_track_20220915_204642.44.csv");
% track_FSM = readmatrix("FSM_CLtrack_20220915_205014.26.csv");

% SL-4 R/B (failed to lock)
% track_TLE = readmatrix("TLE_track_20220915_211033.76.csv");
% track_FSM = readmatrix("FSM_CLtrack_20220915_211604.36.csv");

%% plot TLE
t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

%  az/al vs time
figure(1)
plot(t,track_TLE(:,2:5))
xlabel('Time (sec)')
ylabel('Gimbal Angles (deg)')
legend('measured az','measured el','reference az','reference el')

%  az/al error vs time
% figure(2)
figure;
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
% figure(4)
% plot(t(1:end-1),diff(t))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

% figure(5)
% plot((track_TLE(1:end-1,6:7)/3600.*diff(t)-diff(track_TLE(:,2:3)))*3600)

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

figure;
plot(t_fsm,track_FSM(:,2)*pp,t_fsm,track_FSM(:,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')

figure;
plot(t_fsm,track_FSM(:,4),t_fsm,track_FSM(:,5))
xlabel('Time (sec)')
ylabel('FSM Voltage (V)')

figure;
plot(t_fsm,(track_FSM(:,2)+track_FSM(:,4)/.1)*pp,t_fsm,(track_FSM(:,3)+track_FSM(:,5)/.1)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
% % show time gaps in control step
% figure
% plot(t_fsm(1:end-1),diff(t_fsm))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

figure;
plot(track_FSM(:,1),1e2*track_FSM(:,2:3),track_TLE(:,1),(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)

%% test contol methods
%simulation track (gains 1/0.4, no integral term?)
track_TLE = readmatrix("TLE_track_20220915_204642.44.csv");
t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = track_TLE(1,2:3);%ref_ang(1,:);
% sim_ang(1,:) = ref_ang(1,:)+200/3600;
dt = diff(t);
dv_p = 0;
dv_i = 0;

ref_ang_off = ref_ang; %test target offset
ns = 0; %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);

ref_ang = ref_ang_off;  %test target offset
dt = dt_off;  %test target offset

temp=0;
for i=2:length(ref_rate)

    %reference rate
%     v = 0;
    v = ref_rate(i-1,:);
%     v = (ref_rate(i-1,:)+ref_rate(i,:))/2;

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>0
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
%     k=[1,1];
%     k=[0.4,0.4];
%     k=[0.101,0.101];

    k=[0.4,0.4];

%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));
    
%     temp = (v+1./k.*(dv_p+0*dv_i)); %no lag
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
    temp = (v+1./k.*(dv_p+0*dv_i)); %one step lag

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

%% test
x(1)=1;
x(2)=1;
n=100;
for i=3:n
    tp = -x(i-2)*1/3;

%     if mod(i,2)==0
%         tp = -x(i-2)*0.5;
%     else
%         tp = 0;
%     end

x(i)=x(i-1)+tp;
end
figure()
plot(x)
