function [alpha, delta, lon, lat] = groundTrack_2(t, omega_planet, thetaG0, t0, Y)
% ------------------------------------------------------------------------
% [alpha, delta, lon, lat] = groundTrack_2(t, omega_planet, thetaG0, t0, Y)
% 
% Functin groundTrack_2 computes latitute and longitude, necessary to run the
% function plot_groundTrack starting from the state vector of a orbit. 
%
% INPUTS:
% t              [nx1]       [s]             column vector of time in witch to compute the solution
% omega_planet   [1x1]       [rad/s]         angular velocity of the planet
% thetaG0        [1x1]       [rad]           Greenwich sidderal time at 0 hours UT
% t0             [1x1]       [s]             initial time
% Y              [nx6]       [km, km/s]      state vector: col 1-3 psition vector, col 1-6 velocity vector
% 
% OUTPUTS:
% alpha          [nx1]       [rad]           right ascension
% delta          [nx1]       [rad]           declination
% lon            [nx1]       [deg]           longitude
% lat            [nx1]       [deg]           latitude
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
% Last Update: 24/12/2024
%--------------------------------------------------------------------------

% Separate velocity and position in the state vector
rr = Y(:, 1:3);

% Divide x,y,z components of the position vector
x = rr(:,1);
y = rr(:,2);
z = rr(:,3);

% Calcualte the norm of the position
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