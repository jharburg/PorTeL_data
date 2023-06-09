%% analyze 11/09/22

%% load data
%P_gain 0.4, I_gain 0.2, IR cam: opr 12, centroid 9 pix>250

%CLICK, first light
track_TLE = readmatrix("TLE_track_20221109_214942.09.csv");
track_FSM = readmatrix("FSM_CLtrack_20221109_215220.68.csv");

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
title('CLICK-A Pass Gimbal Angles')
xlim([0 248])

%%%
figure(11)
plot(t,track_TLE(:,2:5))
xlabel('Time (sec)')
ylabel('Gimbal Angles (deg)')
legend('measured az','measured el','reference az','reference el')
title('CLICK-A Pass Gimbal Angles')
xlim([0 248])

figure(12)
plot(t,track_TLE(:,6:7)/3600)
hold on
plot(t,track_TLE(:,8:9))
hold off
xlabel('Time (sec)')
ylabel('Gimbal Rates (deg/sec)')
legend('commanded az rate','commanded el rate','reference az rate','reference el rate')
title('CLICK-A Pass Gimbal Angles Rates')
xlim([0 248])
%%%

%  az/al error (TLE position to mount position) vs time
figure(2)
plot(t,(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
xlim([150 248])
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
title('CLICK-A Pass Gimbal Angles Rates')
xlim([0 248])

% show time gaps in control step
% figure(4)
% plot(t(1:end-1),diff(t))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

tmax = 89;
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
plot([0,100],[-21.75,21.75;-21.75,21.75],'k--')
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
title('CLICK-A Tracking Error')
ylim([0,100])
xlim([0 tmax])
hold on
plot([0,100],[21.75 21.75],'k--')
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
% plot(track_FSM(inds,2)*12.5,track_FSM(inds,3)*12.5,'o')
plot(track_FSM(inds,2)*12.5,track_FSM(inds,3)*12.5,'.')
xlabel('x (um)')
ylabel('y (um)')
title('CLICK-A Laser Centroid Position on PorTeL Tracking Camera')
hold on
x_c = -100:100;
plot(x_c,[sqrt(100^2-x_c.^2);-sqrt(100^2-x_c.^2)],'k--')
hold off
axis equal
legend('Centroid','Requirement')

%%
err = sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp;
1-sum(err>21.75)/length(inds)
% rms(err)