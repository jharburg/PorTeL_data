%star, mount still
track_FSM = readmatrix("FSM_CLtrack_20221018_204445.23.csv");
track_FSM = track_FSM(1:1185,:);

pp=2.7190;  %acrsec/pixel
t_fsm = track_FSM(:,1)-track_FSM(1,1);  %time relative to start of pass

%%
figure;
plot(t_fsm,track_FSM(:,2)*pp,t_fsm,track_FSM(:,3)*pp)
xlabel('Time (sec)')
ylabel('Error (arcsec)')

%%
p1 = polyfit(t_fsm,track_FSM(:,2)*pp,5);
y1 = polyval(p1,t_fsm);
p2 = polyfit(t_fsm,track_FSM(:,3)*pp,5);
y2 = polyval(p2,t_fsm);

figure;
% plot(t_fsm,track_FSM(:,2)*pp,t_fsm,y)
plot(t_fsm,track_FSM(:,2)*pp-y1)
hold on
plot(t_fsm,track_FSM(:,3)*pp-y2)
hold off
xlabel('Time (sec)')
ylabel('Error (arcsec)')

rms(track_FSM(:,2)*pp-y1)
rms(track_FSM(:,3)*pp-y2)