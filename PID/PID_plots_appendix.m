%% Simulate Tracking Algorithm
%% Proportional Control
%load reference positions and rates generated by PorTeL control software
track_TLE = readmatrix("TLE_track_20220628_185120.87.csv");
t = track_TLE(:,1)-track_TLE(1,1);  %set time relative to start of pass

%initialize variables
ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = ref_ang(1,:);
dt = diff(t);
dv_p = 0;
dv_i = 0;

%set time offset
ref_ang_off = ref_ang;
ns = 5;  %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);
ref_ang = ref_ang_off;
dt = dt_off;

%tracking loop
for i=2:(length(ref_rate))

    %reference rate
    v = ref_rate(i-1,:);

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>50
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
    k=[5,5];
    temp = (v+k.*(dv_p+0*dv_i));
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;

end

%plotting
figure(1)
subplot(1,2,1)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional Control')
ylim([-200 200])

%% Proportional-Integral Control

%initialize variables
ref_rate = track_TLE(:,8:9);
ref_ang = track_TLE(:,4:5);
sim_ang = zeros(length(ref_rate),2);
sim_ang(1,:) = ref_ang(1,:);
dt = diff(t);
dv_p = 0;
dv_i = 0;

%set time offset
ref_ang_off = ref_ang;
ns = 5;  %# steps offset (~.2 sec steps)
ref_ang_off(1:end-ns,:) = ref_ang_off(ns+1:end,:);
dt_off = diff(t);
dt_off(1:end-ns,:) = dt_off(ns+1:end,:);
ref_ang = ref_ang_off;
dt = dt_off;

%tracking loop
for i=2:(length(ref_rate))

    %reference rate
    v = ref_rate(i-1,:);

    %correction (proportional term)
    dv_p = (ref_ang(i-1,:)-sim_ang(i-1,:));

    %correction (integral term)
    if max(abs(dv_p))<=100/3600 && i>50
        dv_i = dv_i + dv_p;
    else
        dv_i = 0;
    end

    %update step
    k = [5,5];
    t_i = 0.2;
    temp = (v+k.*(dv_p+0.2/t_i*dv_i));
    sim_ang(i,:) = sim_ang(i-1,:)+dt(i-1)*temp;

end

%plotting
figure(1)
subplot(1,2,2)
plot(t(1:end-ns),(sim_ang(1:end-ns,:)-ref_ang(1:end-ns,:))*3600)
xlabel('Time (sec)')
ylabel('Error (arcsec)')
legend('az','el')
title('Proportional-Integral Control')
ylim([-200 200])
