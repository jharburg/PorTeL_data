%track image analysis
im = imread("trackCam_20221128_193835_967.pgm"); %CLICK centered
im = bitshift(im,-4);

%%
thresh = 1000;
im_t = double(im);
im_t(im_t<thresh) = 0;
y_c = sum(im_t.*(1:256)','all')/sum(double(im_t),'all');
x_c = sum(im_t.*(1:320),'all')/sum(double(im_t),'all');

figure(1)
% imshow(im,[0,thresh],'InitialMagnification','fit')
imshow(im,[],'InitialMagnification','fit')
colorbar
hold on
plot(x_c,y_c,'rx')
hold off

figure(2)
histogram(im,(0:double(max(im(:)))+1)-0.5)

figure(3)
plot(im(round(y_c),:))
hold on
plot(im(:,round(x_c)))
hold off

%% photometry
%constants
c = 2.99792458e8;   %speed of light [m s-2]
h = 6.62607015e-34;  % [J Hz-1]

bg = median(im(:));
n = 7;%APD is size of 16 pixels
im_ob = im(round(y_c)-n:round(y_c)+n,round(x_c)-n:round(x_c)+n)-bg;
flux_c = sum(im_ob(:));  %total counts

figure(4)
% imshow(im,[0,thresh],'InitialMagnification','fit')
imshow(im_ob,[],'InitialMagnification','fit')
colorbar

%opr 9
dt = 4.22e-3;  %exposure time (s)
g = 62;  %nominal gain (e-/count)

flux_p = flux_c*g/dt;  %photons rate (#/s)
P = flux_p*(h*c/1550e-9);  %Power (W)

%% fit gausian

F = @(x,xdata) x(3)/(x(1)*sqrt(2*pi))*exp(-(xdata-x(2)).^2/(2*x(1)^2));

x0_1 = [2 127.5 5e4];
t_1 = ([90: 124 130:160])';
y_1 = double(im(t_1,round(x_c))-bg);

[x_1,resnorm_1,~,exitflag_1,output_1] = lsqcurvefit(F,x0_1,t_1,y_1);

x0_2 = [2 160 5e4];
t_2 = ([135:154 159:175])';
y_2 = double(im(round(y_c),t_2)-bg)';

[x_2,resnorm_2,~,exitflag_2,output_2] = lsqcurvefit(F,x0_2,t_2,y_2);

% x=x0;
figure(5)
t_1p = (1:size(im,1))';
t_2p = (1:size(im,2))';

plot(t_1,y_1,'x',t_1p,F(x_1,t_1p))
hold on
plot(t_2,y_2,'x',t_2p,F(x_2,t_2p))
hold off

% total expected counts
((x_1(3)+x_2(3))/2)*sqrt(2*pi)*((x_1(1)+x_2(1))/2);
