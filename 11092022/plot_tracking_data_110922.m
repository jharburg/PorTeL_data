%% analyze 11/09/22

%% load data
%P_gain 0.4, I_gain 0.2, IR cam: opr 12, centroid 9 pix>250

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

%  az/al error (TLE position to mount position) vs time
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
xlabel('Time (sec)')
ylabel('Gimbal Rates (deg/sec)')
legend('commanded az rate','commanded el rate','reference az rate','reference el rate')

% show time gaps in control step
% figure(4)
% plot(t(1:end-1),diff(t))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

%centroid position
figure(5)
plot(t_fsm,track_FSM(:,2)*pp,t_fsm,track_FSM(:,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
title('CLICK-A Laser Centroid Position on PorTeL Guidance Camera')
legend('x','y')
ylim([-100,100])
hold on
plot([0,100],[-21.75,21.75;-21.75,21.75],'k--')
hold off

%fsm voltage
figure(6)
plot(t_fsm,track_FSM(:,4),t_fsm,track_FSM(:,5))
xlabel('Time (sec)')
ylabel('FSM Voltage (V)')

%centroid position corrected for FSM position
figure(7)
plot(t_fsm,(track_FSM(:,2)+track_FSM(:,4)/.1)*pp,t_fsm,(track_FSM(:,3)+track_FSM(:,5)/.1)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')

% % show time gaps in control step
% figure(8)
% plot(t_fsm(1:end-1),diff(t_fsm))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

figure(9)
plot(track_FSM(:,2)*12.5,track_FSM(:,3)*12.5)
xlabel('x (um)')
ylabel('y (um)')
title('CLICK-A Laser Centroid Position on PorTeL Guidance Camera')
hold on
x_c = -100:100;
plot(x_c,[sqrt(100^2-x_c.^2);-sqrt(100^2-x_c.^2)],'k--')
hold off
axis equal
