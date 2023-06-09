%track image analysis
% im = imread("trackCam_20220513_215734.pgm");
im = imread("trackCam_20220513_215850.pgm");
im = bitshift(im,-4);

%%
figure(1)
imshow(im,[0,4000],'InitialMagnification','fit')
% imshow(im(140:end,140:end),[],'InitialMagnification','fit')

figure(2)
histogram(im,(0:double(max(im(:)))+1)-0.5)
% histogram(im(140:end,140:end),(0:double(max(im(:)))+1)-0.5)

%%
A=double(im(140:end,140:end));
std_A = std(A(:));
mean_A = mean(A(:));
g = mean_A/std_A^2;

%%
figure(3)
plot(im(75,:))