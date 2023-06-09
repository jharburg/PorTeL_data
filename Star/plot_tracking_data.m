%% star tracking for thesis
%% load data
%integral term 0.4, 0.2

%star
track_FSM = readmatrix("FSM_CLtrack_20221018_192544.62.csv");

%% plot FSM
pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

figure;
% plot(track_FSM(:,2:3)*pp)
plot(t_fsm,track_FSM(:,2)*pp,t_fsm,track_FSM(:,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
ylim([-20 20])
xlim([0 130])
hold on
plot([20.5 20.5],[-20 20],'k--')
plot([87.5 87.5],[-20 20],'k--')
hold off
title('Star Tracking Error')
legend('Tracking Cam x-axis','Tracking Cam y-axis')

figure;
plot(t_fsm,sqrt(track_FSM(:,2).^2+track_FSM(:,3).^2)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
ylim([-20 20])
xlim([0 130])
hold on
plot([20.5 20.5],[-20 20],'k--')
plot([87.5 87.5],[-20 20],'k--')
hold off

figure;
plot(t_fsm,track_FSM(:,4),t_fsm,track_FSM(:,5))
xlabel('Time (sec)')
ylabel('FSM Voltage (V)')

% figure;
% plot(t_fsm,(track_FSM(:,2)+track_FSM(:,4)/.1)*pp,t_fsm,(track_FSM(:,3)+track_FSM(:,5)/.1)*pp)
% xlabel('Time (sec)')
% ylabel('Error (arcsec)')

% % show time gaps in control step
% figure
% plot(t_fsm(1:end-1),diff(t_fsm))
% xlabel('Time (sec)')
% ylabel('control step dt (sec)')

% figure;
% plot(track_FSM(:,1),track_FSM(:,2:3),track_TLE(:,1),(track_TLE(:,2:3)-track_TLE(:,4:5))*3600)
% plot(track_FSM(:,1),track_FSM(:,2:3),track_TLE(:,1),500*(track_TLE(:,6:7)/3600-track_TLE(:,8:9)))
%%
i1=find(t_fsm>5,1);
i2=find(t_fsm>20,1);
inds=i1:i2;
rms(track_FSM(inds,2)*pp)
rms(track_FSM(inds,3)*pp)
rms(sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp)

i1=find(t_fsm>22,1);
i2=find(t_fsm>85,1);
inds=i1:i2;
rms(track_FSM(inds,2)*pp)
rms(track_FSM(inds,3)*pp)
rms(sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp)

i1=find(t_fsm>90,1);
i2=find(t_fsm>130,1);
inds=i1:i2;
rms(track_FSM(inds,2)*pp)
rms(track_FSM(inds,3)*pp)
rms(sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp)