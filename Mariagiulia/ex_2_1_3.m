%% ex 2.1.3
% calculate the radius at t = 1.12h
clear all
close all
clc

mu_E = astroConstants(13);
a = 19134.4; %[km]
e = 0.2;
t0 = 0;
theta0 = pi/3;
E0 = 2*atan(sqrt((1-e)/(1+e))*tan(theta0/2)); % 0.8810
%E0 = 0.8810; %[rad]
t = 1.12*60*60;

E = KeplerFunction(t, e, a, mu_E, t0, E0, 1e-6);
theta = 2*atan(sqrt((1+e)/(1-e)*tan(E/2)));
r = a*(1-e^2)/(1+e*cos(theta))

%%
clear all
close all
clc

h = [30738.0, -5367.9, -53520.0]';
om = 11.09; %[deg]
om = om*pi/180;
theta = 306.76;
theta = theta*pi/180;

hnorm = norm(h);
z = [0, 0, 1]';
x = [1, 0, 0]';

N = cross(z,h)/norm(cross(z,h))
OMEGA = acos(dot(N,x))*180/pi
i = acos(h(3)/hnorm)*180/pi

%% 
r = 1e+3*[-3.8122; -1.6241; -5.4898];
normr = norm(r)
rxy = norm([r(1) r(2)])