i1=find(t_fsm>2,1);
i2=find(t_fsm>111,1);
inds=i1:i2;
pp=2.7190;  %acrsec/pixel
% rms(track_FSM(inds,2)*pp)
% rms(track_FSM(inds,3)*pp)
rms(sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp)

% rms((track_FSM(inds,2)+track_FSM(inds,4)/.1)*pp)
% rms((track_FSM(inds,3)+track_FSM(inds,5)/.1)*pp)
rms(sqrt((track_FSM(inds,2)+track_FSM(inds,4)/.1).^2 ...
    +(track_FSM(inds,3)+track_FSM(inds,5)/.1).^2)*pp)

% figure()
% histogram(track_FSM(inds,2)*pp)

err=sqrt(track_FSM(inds,2).^2+track_FSM(inds,3).^2)*pp;
1-sum(err>21.75)/length(inds)

err=sqrt((track_FSM(inds,2)+track_FSM(inds,4)/.1).^2 ...
    +(track_FSM(inds,3)+track_FSM(inds,5)/.1).^2)*pp;
1-sum(err>21.75)/length(inds)