%track image analysis
im = imread("trackCam_20221109_215300_257.pgm");
% im = imread("trackCam_20221109_215216_092.pgm");
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
n = 20;
im_ob = im(round(y_c)-n:round(y_c)+n,round(x_c)-n:round(x_c)+n)-bg;
flux_c = sum(im_ob(:));  %total counts

figure(4)
% imshow(im,[0,thresh],'InitialMagnification','fit')
imshow(im_ob,[],'InitialMagnification','fit')
colorbar

%opr 12
dt = 11.67e-3;  %exposure time (s)
g = 11;  %nominal gain (e-/count)

flux_p = flux_c*g/dt;  %photons rate (#/s)
P = flux_p*(h*c/1550e-9);  %Power (W)
