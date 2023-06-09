%% test contol methods
%simulation track (gains 1/0.4, no integral term)
track_TLE = readmatrix("TLE_track_20220628_185120.87.csv");

t = track_TLE(:,1)-track_TLE(1,1);  %time relative to start of pass

ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = ref_ang(1,:);
% sim_ang(1,:) = track_TLE(1,2:3);%ref_ang(1,:);
% sim_ang(1,:) = ref_ang(1,:)+200/3600;
dt = diff(t);
dv_p = 0;
dv_i = 0;

ref_ang_off = ref_ang; %test target offset
ns = 5; %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);

ref_ang = ref_ang_off;  %test target offset
dt = dt_off;  %test target offset

temp=0;
for i=2:(length(ref_rate))

    %reference rate
%     v = 0;
    v = ref_rate(i-1,:);
%     v = (ref_rate(i-1,:)+ref_rate(i,:))/2;

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>50
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
%     k=[1,1];
%     k=[0.4,0.4];
%     k=[0.2,0.2];
%     k=[0.105,0.105];

    k=[.2,.2];

%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));

    temp = (v+1./k.*(dv_p+0*dv_i)); %no lag
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
%     temp = (v+1./k.*(dv_p+0.25*dv_i)); %one step lag

end

% %  az/al vs time
% figure()
% plot(t,[sim_ang ref_ang])
% xlabel('Time (sec)')
% ylabel('Gimbal Angles (deg)')
% legend('sim az','sim el','reference az','reference el')

figure(1)
subplot(1,2,1)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional Control')
ylim([-200 200])

%%%
figure(2)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional Control')
ylim([-200 200])
%%%

%test
er1=sim_ang-ref_ang;
% rms(er1(1000:1700,:)*3600)
rms(er1(100:end-100,:)*3600)
% rms(er1*3600)

max(abs(er1(100:end-100,:)*3600))

%%
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
ns = 5; %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);

ref_ang = ref_ang_off;  %test target offset
dt = dt_off;  %test target offset

temp=0;
for i=2:(length(ref_rate))

    %reference rate
%     v = 0;
    v = ref_rate(i-1,:);
%     v = (ref_rate(i-1,:)+ref_rate(i,:))/2;

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>50
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
%     k=[1,1];
%     k=[0.4,0.4];
%     k=[0.2,0.2];
%     k=[0.105,0.105];

    k=[.2,.2];

%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));

    temp = (v+1./k.*(dv_p+1*dv_i)); %no lag
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
%     temp = (v+1./k.*(dv_p+0.25*dv_i)); %one step lag

end

figure(1)
subplot(1,2,2)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional-Integral Control')
ylim([-200 200])

%%%
figure(3)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional-Integral Control')
ylim([-200 200])
%%%

%test
er2=sim_ang-ref_ang;
% rms(er2(1000:1700,:)*3600)
rms(er2(100:end-100,:)*3600)
% rms(er2*3600)

max(abs(er2(100:end-100,:)*3600))

%%
%% sim over t_i's
t_is = 0.1:0.1:1;
max_err = zeros(2,length(t_is));
rms_err = zeros(2,length(t_is));

for j=1:length(t_is)
t_i=t_is(j);
ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = track_TLE(1,2:3);%ref_ang(1,:);
% sim_ang(1,:) = ref_ang(1,:)+200/3600;
dt = diff(t);
dv_p = 0;
dv_i = 0;

ref_ang_off = ref_ang; %test target offset
ns = 5; %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);

ref_ang = ref_ang_off;  %test target offset
dt = dt_off;  %test target offset

temp=0;
for i=2:(length(ref_rate))

    %reference rate
%     v = 0;
    v = ref_rate(i-1,:);
%     v = (ref_rate(i-1,:)+ref_rate(i,:))/2;

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>50
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
%     k=[1,1];
%     k=[0.4,0.4];
%     k=[0.2,0.2];
%     k=[0.105,0.105];

    k=[5,5];

%     sim_ang(i,:) = sim_ang(i-1,:)+0.2028*(v+1./k.*(dv_p+1*dv_i));

    temp = (v+k.*(dv_p+0.2./t_i*dv_i)); %no lag
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;
%     temp = (v+1./k.*(dv_p+0.25*dv_i)); %one step lag

end

% figure()
% plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
% xlabel('Time (sec)')
% ylabel('Error (arcsec)')
% legend('az','el')
% title('Proportional-Integral Control')
% ylim([-200 200])

%test
er2=sim_ang-ref_ang;
% rms(er2(1000:1700,:)*3600)
rms_err(:,j)=rms(er2(100:end-100,:)*3600);
% rms(er2*3600)

max_err(:,j)=max(abs(er2(100:end-100,:)*3600));

end
figure
plot(t_is,max_err,'o')
figure
plot(t_is,rms_err,'o')
