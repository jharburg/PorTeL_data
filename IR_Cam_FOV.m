% telescope focused on IR cam (no extra optics)
p1=rad2deg(atan(1/2*12.5e-6/2.8))*3600*2;  %arcsec/pixel 

% telescope, collimator, lens, IR cam
th2 = atan(1/2*12.5e-6/.0254); %theta 2 (rad)
x2 = tan(th2)*.075;  % (m)
p2 = rad2deg(atan(x2/2.8))*3600*2;  %%arcsec/pixel

% telescope, collimator, lens, APD
th2 = atan(1/2*200e-6/.0254); %theta 2 (rad)
x2 = tan(th2)*.075;  % (m)
apd = rad2deg(atan(x2/2.8))*3600*2;  %%arcsec/pixel

%% FSM walk
f1 = 2.8;  %scope focal length (m)
f2 = 0.075;  %collimating lens focal length (m)
d = 0.1;  %dist from collimating lens to fsm (m)

th = deg2rad(0:.01:.1);  %angle of arrival (rad)

h_1 = f1*tan(th);  %off axis radius at scope focus
h_2 = (f1+f2)*tan(th)-d*tan(th*f1/f2); %off axis radius at d

figure()
plot(rad2deg(th),[h_1;h_2])

ap = 0.28;  %aperture diameter (m)
d2 = 0.1;  %dist focus to fsm (original setup, m)
h_1p = (f1-d2)*tan(th)+ap*d2/f1/2; %off axis max at fsm
h_2p = (f1+f2)*tan(th)-d*tan(th*f1/f2)+ap*f2/f1/2; %off axis max at fsm
h_2m = (f1+f2)*tan(th)-d*tan(th*f1/f2)-ap*f2/f1/2; %off axis min at fsm

figure()
plot(rad2deg(th),[h_1p;h_2p;h_2m])

%% tel plot old
f1 = 2.8;  %scope focal length (m)
ap1 = .28;  %scope aperture diameter
fsm_d = 0.026;  %fsm diameter
ir_l = 0.1;  %dist fsm-ir  (IR cam at focus)

figure(10)
%scope body
plot([0;f1],[ap1/2 -ap1/2;ap1/2 -ap1/2],'k')
hold on
axis equal
%fsm
plot([f1-ir_l;f1-ir_l],[fsm_d/2;-fsm_d/2],'k')

%in rays
x=(0:.005:f1);
th=deg2rad(0.15);  %incoming angle (rad)

%r1, top edge aperture
r1x = x;
r1y = ap1/2+(f1*tan(th)-ap1/2)/f1*x;
r1x = [-.5 0 r1x];
r1y = [ap1/2-tan(th)*.5 ap1/2 r1y];
plot(r1x,r1y,'b')
%r2, bottom edge aperture
r2x = x;
r2y = -ap1/2+(f1*tan(th)+ap1/2)/f1*x;
r2x = [-.5 0 r2x];
r2y = [-ap1/2-tan(th)*.5 -ap1/2 r2y];
plot(r2x,r2y,'b')
%r3, center aperture
r3x = x;
r3y = (f1*tan(th))/f1*x;
r3x = [-.5 0 r3x];
r3y = [-tan(th)*.5 0 r3y];
plot(r3x,r3y,'b')

%
hold off

%% tel plot new
f1 = 2.8;  %scope focal length (m)
f2 = 0.075;  %collimating lens focal length (m)
ap1 = .28;  %scope aperture diameter
f2_d = 0.026;  %relay lens diameter
fsm_d = 0.026;  %fsm diameter
ir_d = 0.026;  %IR cam diameter
fsm_l = 0.02;  %dist f2-fsm 
ir_l = 0.08;  %dist fsm-ir 

figure(11)
%scope body
plot([0;f1],[ap1/2 -ap1/2;ap1/2 -ap1/2],'k')
hold on
axis equal
%relay
plot([f1+f2;f1+f2],[f2_d/2;-f2_d/2],'k')
%fsm
plot([f1+f2+fsm_l;f1+f2+fsm_l],[fsm_d/2;-fsm_d/2],'k')
%IR
plot([f1+f2+fsm_l+ir_l;f1+f2+fsm_l+ir_l],[ir_d/2;-ir_d/2],'k')


%in rays
x=(0:.005:f1+f2);
x2=(f1+f2:0.005:f1+f2+fsm_l+ir_l);
% th=deg2rad(0);  %incoming angle (rad)
th=deg2rad([0.5;0;-0.5]);  %incoming angle (rad)
colororder({'b','r','g'})

%r1, top edge aperture
r1x = x;
r1y = ap1/2+(f1*tan(th)-ap1/2)/f1.*x;
r1x = [-.5 0 r1x];
r1y = [ap1/2-tan(th)*.5 repmat(ap1/2,size(th)) r1y];
r1x = [r1x x2];
r1y = [r1y r1y(:,end)-tan(th*f1/f2).*(x2-x2(1))];
plot(r1x,r1y)
%r2, bottom edge aperture
r2x = x;
r2y = -ap1/2+(f1*tan(th)+ap1/2)/f1.*x;
r2x = [-.5 0 r2x];
r2y = [-ap1/2-tan(th)*.5 repmat(-ap1/2,size(th)) r2y];
r2x = [r2x x2];
r2y = [r2y r2y(:,end)-tan(th*f1/f2).*(x2-x2(1))];
plot(r2x,r2y)
%r3, center aperture
r3x = x;
r3y = (f1*tan(th))/f1.*x;
r3x = [-.5 0 r3x];
r3y = [-tan(th)*.5 repmat(0,size(th)) r3y];
r3x = [r3x x2];
r3y = [r3y r3y(:,end)-tan(th*f1/f2).*(x2-x2(1))];
plot(r3x,r3y)

%
hold off

