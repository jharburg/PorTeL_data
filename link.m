%% optic link
%% contants
kb = 1.380649e-23;   %Boltzmann constant [J K-1]
c = 2.99792458e8;   %speed of light [m s-2]
h = 6.62607015e-34;  % [J Hz-1]

R_s = 6.9634e8;  %sun radius [m]
au = 1.496e11;  %1au [m]

%% BB
T = 5800; %temp [k]
l = (1:1:5000)*1e-9; %wavelength [m]
B = 2*h*c^2./l.^5*1./(exp(h*c./(l*kb*T))-1);

figure()
plot(l,B)

trapz(l,B)*pi*R_s^2/au^2;

%% laser
l = 1550e-9;  %laser wavelength (m)
P_tx = 200e-3;  %tx power (W)
pnt_er = 0;  %pointing error (rad)
% dist = 500e3;  %dist to target (m)
dist = 517.6e3;
r_rx = .28/2;  %recieve aperture radius (m)

FWHM_div_tx= 1300e-6;  %FWHM beam divergence(rad)
rad_tx = l/(pi*(FWHM_div_tx/sqrt(2*log(2))));  % beam waist w_0 (m)

rad_rx = rad_tx*sqrt(1+(l*dist/(pi*rad_tx^2))^2);  %beam radius w at dist (m)
dist_off_rx = tan(pnt_er)*dist;  %off axis dist (m)
irad_rx = P_tx*2/(pi*rad_rx^2)*exp(-2*dist_off_rx^2/rad_rx^2);  %rx iradiance (W/m^2)
P_rx = irad_rx*pi*r_rx^2;  %rx power (W)

l_other = 11.27+3.62;  %other losses
P_rx = P_rx/(10^(l_other/10));

n_phot_l = P_rx/(h*c/l);  %photon rate (#/s)

%% Vega
Flux_V = 3.06e-10;  %vega J band flux (erg s-1 cm-2 A-1)
Flux_V = Flux_V*10^-7*100^2*10^10;  %(W m-3)
bw = 1700e-9 - 900e-9;  %SWIR cam bandwidth, no filter (m)
% bw = 10e-9;  %SWIR cam bandwidth,  filter (m)
irad_rx = Flux_V*bw;  %rx iradiance (W/m^2)
P_rx = irad_rx*pi*r_rx^2;  %rx power (W)

n_phot_v = P_rx/(h*c/1550e-9);  %photon rate (#/s)

%% Solar-type
T = 5800; %temp [k]
l = (900:1:1700)*1e-9; %wavelength, no filter [m]
% l = (1550-6:1:1550+6)*1e-9; %wavelength, filter [m]
B = 2*h*c^2./l.^5*1./(exp(h*c./(l*kb*T))-1);  %spectral radiance (W/m^3/sr)
B_phot = B./(h*c./l);  %photon flux (#/s/m^3/sr)

% figure()
% plot(l,B_phot)

B_phot = B_phot*pi*R_s^2/au^2;  %photon flux (#/s/m^3)
B_phot_0 = B_phot*10^(-0.4*26.74);  %photon flux of mag 0 solar-type star (#/s/m^3)
n_phot_0 = trapz(l,B_phot_0)*pi*r_rx^2;  %photon rate (#/s)

B_0 = trapz(l,B_phot_0.*(h*c./l))*pi*r_rx^2;  %Power through aperature (W)

V_laser = -2.5*log10(n_phot_l/n_phot_0);

%% sky background
sky = 10^6;  %sky background est (#/s/m^2/um/arcsec^2)
sky = sky*10^6;  %sky background est (#/s/m^3/arcsec^2)
sky = sky*2.719^2*pi*r_rx^2;  %sky background est (#/s/m/pix)

% bw = 1700e-9 - 900e-9;  %SWIR cam bandwidth, no filter (m)
bw = 10e-9;  %SWIR cam bandwidth,  filter (m)

n_phot_b = bw*sky;  %background photon rate (#/s/pix)

%% NEI
nei = 1.2e9;  %(#/cm^2/s)
nei = nei*(12.5e-4)^2;  %(#/pix/s)

%% ir cam im sim
% n_phot_l

fwhm = 2;  %in pixels
sigma = fwhm/(sqrt(8*log(2)));
A = 1;  %total signal
% figure()
% inds = -10:0.1:10;
% plot(inds,normpdf(inds,0,sigma))

peak = A/(2*pi*sigma^2);  %peak pixel center val
peak_pix = A*erf(1/(2*sqrt(2)*sigma))^2;  %peak pixel total val (analytical), source at center
% peak_pix_p = A/4*erf(1/(sqrt(2)*sigma))^2;  %peak pixel total val (analytical), source at corner

%arbitrary limits
% xmin = -.5;xmax=0.5;
% ymin = -.5;ymax=0.5;
% peak_pix = A/4*(erf(xmax/(sqrt(2)*sigma))-erf(xmin/(sqrt(2)*sigma)))*...
%     (erf(ymax/(sqrt(2)*sigma))-erf(ymin/(sqrt(2)*sigma)));

%% oprs
opr_dts = 1e-3*[0.20 0.28 0.39 0.55 0.78 1.09 1.53 2.14 3.01 4.22 5.92 8.31 11.67 16.37];  %(s)
opr_gs = [11000 6200 3500 2000 1100 620 350 190 110 62 35 19 11 6.1];  %(e-/count)

opr_n = 6;
opr_dt = opr_dts(opr_n+1);
opr_g = opr_gs(opr_n+1);

n_phot = n_phot_l;
% n_phot = n_phot_0*0.1;

peak_pix_adu = n_phot*opr_dt/opr_g*peak_pix;

offset = 2;  %pix
offset_pix_adu = n_phot*opr_dt/opr_g*(exp(-offset^2/(2*sigma^2))/(2*pi*sigma^2));  %val at center of pix



