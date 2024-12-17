%% LAB chapeter 2
% 30 october
%
% ex 1: computation of ground tracks

%% 1
clear all
close all
clc

a = 8350; 
e = 0.1976;
i = 60/180 * pi; %[deg->rad]
OMEGA = 270/180 * pi; %[deg->rad]
omega = 45/180 * pi; %[deg->rad]
theta0 = 230/180 * pi; %[deg->rad]
r0 = [-4578.219; -801.084; -7929.708];
v0 = [0.800; -6.037; 1.385];
thetaG0 = 0;
t0 = 0;

mu_E = astroConstants(13);
T = 2*pi*sqrt(a^3/mu_E);
omegaE = 15.04/180*pi/60/60; 
t = linspace(t0, 3.25*T, 10000);
t = t';

%[rr,vv] = par2car(a,e,i,OMEGA,omega,theta0, mu_E)
%[aa, ee, ii, OO, oo, thth] = car2par(r0, v0, mu_E)
[alpha, delta, lon, lat] = groundTrack('Car', t, omegaE, thetaG0, mu_E, t0, [r0; v0]);

%[alpha, delta, lon, lat] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a,e,i,OMEGA,omega,theta0]);

lat = lat./pi .* 180;
lon = lon./pi .* 180;

plotGroundTrack(lon, lat, t)

%%
clear all
close all
clc

a = 26600; 
e = 0.74;
i = 63.4/180 * pi; %[deg->rad]
OMEGA = 50/180 * pi; %[deg->rad]
omega = 280/180 * pi; %[deg->rad]
theta0 = 0/180 * pi; %[deg->rad]
%r0 = [-4578.219; -801.084; -7929.708];
%v0 = [0.800; -6.037; 1.385];
thetaG0 = 0;
t0 = 0;

mu_E = astroConstants(13);
T = 2*pi*sqrt(a^3/mu_E);
omegaE = 15.04/180*pi/60/60; 
t = linspace(t0, 30*T, 10000);
t = t';

[alpha, delta, lon, lat] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a,e,i,OMEGA,omega,theta0]);

lat = lat./pi .* 180;
lon = lon./pi .* 180;

plotGroundTrack(lon, lat, t)

%%
clear all
close all
clc

k = 12;
m = 1;
a_gt = a_groundTrack(k, m);
a = 8350; 
e = 0.1976;
i = 60/180 * pi; %[deg->rad]
OMEGA = 270/180 * pi; %[deg->rad]
omega = 45/180 * pi; %[deg->rad]
theta0 = 230/180 * pi; %[deg->rad]
r0 = [-4578.219; -801.084; -7929.708];
v0 = [0.800; -6.037; 1.385];
thetaG0 = 0;
t0 = 0;

mu_E = astroConstants(13);
T = 2*pi*sqrt(a_gt^3/mu_E);
omegaE = 15.04/180*pi/60/60; 
t = linspace(t0, T, 10000);
t = t';


[alpha, delta, lon_gt, lat_gt] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a_gt,e,i,OMEGA,omega,theta0]);
[alpha, delta, lon, lat] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a,e,i,OMEGA,omega,theta0]);

lat_gt = lat_gt./pi .* 180;
lon_gt = lon_gt./pi .* 180;
lat = lat./pi .* 180;
lon = lon./pi .* 180;


% Carica l'immagine
img = imread('EarthTexture.jpg');

% Ottieni le dimensioni dell'immagine
%[rows, cols, ~] = size(img);

% Crea il plot
figure;
ax = axes;

% Mostra l'immagine come sfondo, ridimensionata per coprire l'intero range degli assi
imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', ax);
set(ax, 'YDir', 'normal');
hold on;
axis on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

% Imposta i limiti degli assi per coprire l'immagine
xlim([-180 180]);
ylim([-90 90]);

% Tracciare una griglia per i meridiani e paralleli (opzionale)
plot([0 0], ylim, 'k', 'LineWidth', 1); % Meridiano di Greenwich
plot(xlim, [0 0], 'k', 'LineWidth', 1); % Equatore

plot(lon_gt, lat_gt, '*b');
plot(lon_gt(1), lat_gt(1), 'or', LineWidth=2);
plot(lon_gt(end), lat_gt(end), 'ob', LineWidth=2);
plot(lon, lat, '*r');
plot(lon(1), lat(1), 'or', LineWidth=2);
plot(lon(end), lat(end), 'ob', LineWidth=2);

%%
clear all
close all
clc

k = 2;
m = 1;
a_gt = a_groundTrack(k, m);
a = 26600; 
e = 0.74;
i = 63.4/180 * pi; %[deg->rad]
OMEGA = 50/180 * pi; %[deg->rad]
omega = 280/180 * pi; %[deg->rad]
theta0 = 0/180 * pi; %[deg->rad]
%r0 = [-4578.219; -801.084; -7929.708];
%v0 = [0.800; -6.037; 1.385];
thetaG0 = 0;
t0 = 0;

mu_E = astroConstants(13);
T = 2*pi*sqrt(a_gt^3/mu_E);
omegaE = 15.04/180*pi/60/60; 
t = linspace(t0, 100*T, 10000);
t = t';


[alpha, delta, lon_gt, lat_gt] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a_gt,e,i,OMEGA,omega,theta0]);
[alpha, delta, lon, lat] = groundTrack('Kep', t, omegaE, thetaG0, mu_E, t0, [a,e,i,OMEGA,omega,theta0]);

lat_gt = lat_gt./pi .* 180;
lon_gt = lon_gt./pi .* 180;
lat = lat./pi .* 180;
lon = lon./pi .* 180;


% Carica l'immagine
img = imread('EarthTexture.jpg');

% Ottieni le dimensioni dell'immagine
%[rows, cols, ~] = size(img);

% Crea il plot
figure;
ax = axes;

% Mostra l'immagine come sfondo, ridimensionata per coprire l'intero range degli assi
imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', ax);
set(ax, 'YDir', 'normal');
hold on;
axis on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

% Imposta i limiti degli assi per coprire l'immagine
xlim([-180 180]);
ylim([-90 90]);

% Tracciare una griglia per i meridiani e paralleli (opzionale)
plot([0 0], ylim, 'k', 'LineWidth', 1); % Meridiano di Greenwich
plot(xlim, [0 0], 'k', 'LineWidth', 1); % Equatore

plot(lon_gt, lat_gt, 'ob');
plot(lon_gt(1), lat_gt(1), 'or', LineWidth=2);
plot(lon_gt(end), lat_gt(end), 'ob', LineWidth=2);
plot(lon, lat, '*r');
plot(lon(1), lat(1), 'or', LineWidth=2);
plot(lon(end), lat(end), 'ob', LineWidth=2);






