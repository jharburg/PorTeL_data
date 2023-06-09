%% analyze 11/09/22

%% load data
%testing integral term 0.4, 0.2, opr 12

%CLICK, first light
track_TLE = readmatrix("TLE_track_20221109_214942.09.csv");
track_FSM = readmatrix("FSM_CLtrack_20221109_215220.68.csv");
%% plot TLE
t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass
% t = datetime(track_TLE(:,1),'ConvertFrom','posixtime');

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

% figure()
% plot(t,track_TLE(:,6:7)/3600-track_TLE(:,8:9))
% hold on
% plot(t,track_TLE(:,10:11)/100)
% hold off

figure()
%10 is elevation axis
% plot(t,track_TLE(:,10:11),t,track_TLE(:,12:13))
plot(t,track_TLE(:,10)-(track_TLE(:,12)-160))
% plot(t,track_TLE(:,10),t,(track_TLE(:,12)-160))

% show time gaps in control step
% figure(4)
% plot(t(1:end-1),diff(t))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

figure(5)
plot(t(2:end),(track_TLE(1:end-1,6:7)/3600.*diff(t)-diff(track_TLE(:,2:3)))*3600)
% plot(t(1:end-1),track_TLE(1:end-1,6:7).*diff(t),t(1:end-1),diff(track_TLE(:,2:3))*3600)

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

figure;
plot(t_fsm,track_FSM(:,2)*pp,t_fsm,track_FSM(:,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
title('CLICK-A Laser Centroid Position on PorTeL Guidance Camera')
legend('x','y')
ylim([-100,100])
hold on
plot([0,100],[-21.75,21.75;-21.75,21.75],'k--')
hold off

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
% plot(track_FSM(:,1),track_FSM(:,2:3),track_TLE(:,1),(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)
plot(track_FSM(:,1),track_FSM(:,2:3),track_TLE(:,1),500*(track_TLE(:,6:7)/3600-track_TLE(:,8:9)))

figure;
plot(track_FSM(:,2)*12.5,track_FSM(:,3)*12.5)
xlabel('x (um)')
ylabel('y (um)')
title('CLICK-A Laser Centroid Position on PorTeL Guidance Camera')
hold on
x_c = -100:100;
plot(x_c,[sqrt(100^2-x_c.^2);-sqrt(100^2-x_c.^2)],'k--')
hold off
axis equal

%% analyse oscillation
figure()
plot(t,(track_TLE(:,7)-track_TLE(:,9)*3600)*0.2)%dist comanded to move in 1 step
hold on
plot(t,track_TLE(:,10)*pp)%centroidw/o fsm
plot(t,(track_TLE(:,12)-160)*pp)%centroid
plot(t(1:end-1),-diff(track_TLE(:,10)*pp))
plot(t(1:end-1),-diff((track_TLE(:,12)-160)*pp))
hold off


%% test contol methods
%simulation track
track_TLE = readmatrix("TLE_track_20221109_214942.09.csv");

t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
% sim_ang(1,:) = track_TLE(1,2:3);%ref_ang(1,:);
sim_ang(1,:) = ref_ang(1,:)+200/3600;
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
    if max(abs(dv_p))<=100/3600 && i>20
        if 1%mod(i,3)==0
            dv_i = dv_i + dv_p;
        end
    else
        dv_i = 0;
    end

    %update step
%     k=[1,1];
%     k=[0.4,0.4];
%     k=[0.101,0.101];

    k=[0.4,0.4];

%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));

%     temp = (v+1./k.*(dv_p+0.2*dv_i)); %no lag

    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
    if abs(sim_ang(i,1)-ref_ang(i,1))>180
        sim_ang(i,1) = sim_ang(i,1)-360*sign(sim_ang(i,1)-ref_ang(i,1));
    end

    temp = (v+1./k.*(dv_p+0.2*dv_i)); %one step lag

end

% %  az/al vs time
% figure()
% plot(t,[sim_ang ref_ang])
% xlabel('Time (sec)')
% ylabel('Gimbal Angles (deg)')
% legend('sim az','sim el','reference az','reference el')

figure()
err = sim_ang-ref_ang;
% err(:,1) = err(:,1).*cos(deg2rad(sim_ang(:,2)));%plots xy on track cam err
plot(t,err*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')

% figure()
% plot(t,(ref_ang-track_TLE(:,4:5))*3600)

