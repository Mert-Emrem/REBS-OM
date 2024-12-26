function [alpha, delta, lon, lat] = groundTrack_2(t, omega_planet, thetaG0, t0, Y)
% Functin groundTrack_2 computes latitute and longitude, necessary to run the
% function plot_groundTrack starting from the state vector of a orbit. 
%
% STRUCTURE:
%   [alpha, delta, lon, lat] = groundTrack_2(t, omega_planet, thetaG0, t0, Y)
%
% INPUT:
%  - t              [nx1]   column vector of time in witch to compute the solution
%  - omega_planet   [1x1]   angular velocity of the planet
%  - thetaG0        [1x1]   Greenwich sidderal time at 0 hours UT
%  - t0             [1x1]   initial time
%  - Y              [nx6]   state vector
%
% OUTPUT:
%  - alpha          [nx1]   right ascension
%  - delta          [nx1]   declination
%  - lon            [nx1]   longitude
%  - lat            [nx1]   latitude
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
%--------------------------------------------------------------------------

% Separate velocity and position in the state vector
rr = Y(:, 1:3);
%vv = Y(:, 4:6);

% Divide x,y,z components of the position vector
x = rr(:,1);
y = rr(:,2);
z = rr(:,3);

% Calcualte the norm ot the position
rnorm = vecnorm(rr, 2, 2);

% Calculation of alpha delta longitude and latitude
alpha = atan2(y,x);
delta = asin(z./rnorm);
thetaG = thetaG0 + omega_planet*(t - t0);

lon = alpha - thetaG;

for i = 1:length(lon)
    while lon(i)>pi
        lon(i) = lon(i) - 2*pi;
    end
    while lon(i)<-pi
        lon(i) = lon(i) + 2*pi;
    end
end

lat = delta;

% Convert latitude and longitude from radiants to degrees
lat = lat./pi .* 180;
lon = lon./pi .* 180;

end