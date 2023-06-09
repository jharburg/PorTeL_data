%% analyze 11/19/22

%% load data
% integral term 0.4, 0.2

%ISS
track_TLE = readmatrix("TLE_track_20221119_180357.47.csv");
track_FSM = readmatrix("FSM_CLtrack_20221119_180732.23.csv");

track_TLE(:,2)=mod(track_TLE(:,2),360);
track_TLE(:,4)=mod(track_TLE(:,4),360);

%% plot TLE
t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass
% t = datetime(track_TLE(:,1),'ConvertFrom','posixtime');

%  az/al vs time
figure(1)
subplot(1,2,1)
plot(t,track_TLE(:,2:5))
xlabel('Time (sec)')
ylabel('Gimbal Angles (deg)')
legend('measured az','measured el','reference az','reference el')
title('ISS Pass Gimbal Angles')
% xlim([0 248])

%  az/al error (TLE position to mount position) vs time
figure(2)
plot(t,(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
% xlim([150 248])
title('Difference in Mount Position from TLE estimated Position')

%  az/al rate vs time
% figure(3)
figure(1)
subplot(1,2,2)
plot(t,track_TLE(:,6:7)/3600)
hold on
plot(t,track_TLE(:,8:9))
hold off
xlabel('Time (sec)')
ylabel('Gimbal Rates (deg/sec)')
legend('commanded az rate','commanded el rate','reference az rate','reference el rate')
title('ISS Pass Gimbal Angles Rates')
% xlim([0 248])

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

tmax = 111;
i1=find(t_fsm>=0,1);
i2=find(t_fsm>tmax,1);
inds=i1:i2;

%centroid position
figure(5)
plot(t_fsm(inds),track_FSM(inds,2)*pp,t_fsm(inds),track_FSM(inds,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
title('CLICK-A Tracking Error')
ylim([-100,100])
xlim([0 tmax])
hold on
plot([0,120],[-21.75,21.75;-21.75,21.75],'k--')
hold off
legend('tracking cam x-axis','tracking cam y-axis','Requirement')
% colororder({'k'})
% yyaxis right
% ylim([-100,100]/pp*12.5)
% ylabel('Error (um)')

figure(6)
plot(t_fsm(inds),sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
title('ISS Tracking Error')
ylim([0,100])
xlim([0 tmax])
hold on
plot([0,120],[21.75 21.75],'k--')
hold off
legend('Error','Requirement')
colororder({'k'})
yyaxis right
ylim([0,100]/pp*12.5)
ylabel('Error (um)')

%fsm voltage
figure(7)
plot(t_fsm,track_FSM(:,4),t_fsm,track_FSM(:,5))
xlabel('Time (sec)')
ylabel('FSM Voltage (V)')

%centroid position corrected for FSM position
figure(8)
plot(t_fsm,(track_FSM(:,2)+track_FSM(:,4)/.1)*pp,t_fsm,(track_FSM(:,3)+track_FSM(:,5)/.1)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')

% % show time gaps in control step
% figure(9)
% plot(t_fsm(1:end-1),diff(t_fsm))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

figure(10)
plot(track_FSM(inds,2)*12.5,track_FSM(inds,3)*12.5,'o')
xlabel('x (um)')
ylabel('y (um)')
title('ISS Centroid Position on PorTeL Tracking Camera')
hold on
x_c = -100:100;
plot(x_c,[sqrt(100^2-x_c.^2);-sqrt(100^2-x_c.^2)],'k--')
hold off
axis equal
legend('Centroid','Requirement')

%% test contol methods
%simulation track
track_TLE = readmatrix("TLE_track_20221018_195247.56.csv");

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

    temp = (v+1./k.*(dv_p+0*dv_i)); %no lag
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
%     temp = (v+1./k.*(dv_p+0*dv_i)); %one step lag

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
