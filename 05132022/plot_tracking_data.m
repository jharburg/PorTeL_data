%% analyze 05/13/22 ISS pass tracking data

%% load data
track_TLE = readmatrix("TLE_track_20220513_215206.30.csv");
track_fsm = readmatrix("FSM_CLtrack_20220513_215840.32.csv");

%% plot TLE
t=track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

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
xlabel('Time (sec)')
ylabel('Gimbal Rates (deg/sec)')
legend('commanded az rate','commanded el rate','reference az rate','reference el rate')

% show time gaps in control step
figure(4)
plot(t(1:end-1),diff(t))
xlabel('Time (sec)')
ylabel('control step dt (sec)')

%% plot fms
t_fsm = track_fsm(:,1)-track_fsm(1,1);  %time relative to start of pass

% plot fsm
figure(5)
p2=2.7190;  %IR cam arcsec/pix
plot(t_fsm,track_fsm(:,2:3)*p2)
xlabel('Time (sec)')
ylabel('centroid position on guidance tracker (arcsec)')
legend('x','y')

figure(6)
plot(track_fsm(:,2)*p2,track_fsm(:,3)*p2,'x')
xlabel('centroid x position on guidance tracker  (arcsec)')
ylabel('centroid y position on guidance tracker  (arcsec)')
axis equal

figure(7)
plot(t_fsm,track_fsm(:,4:5))
xlabel('Time (sec)')
ylabel('FSM Voltage (V)')
legend('x','y')

% figure(8)
% plot(t_fsm(1:end-1),diff(t_fsm))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')