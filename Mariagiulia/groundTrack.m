function [alpha, delta, lon, lat] = groundTrack(type, t, omegaE, thetaG0, mu, t0, coord)

% ground track function to calulate latitute and longitude
%
% INPUT:
% type      [str] cartesian data or keplerian data  - Kep or Car
% t         [nx1] column vector of time in witch to compute the solution
% omegaE    [1x1] angular velocity of the Earth
% thetaG0   [1x1] Greenwich sidderal time at 0 hours UT
% mu        [1x1] gravitational constant
% t0        [1x1] initial time
% coord     
% - if Kep: [6x1] with a, e, i, OMEGA, omega, theta
% - if Car: [1x6] with r0 (1,2,3) and v0 (4,5,6)
%
% OUTPUT
% alpha     [nx1] right ascension
% delta     [nx1] declination
% lon       [nx1] longitude
% lat       [nx1] latitude

if strcmp(type, 'Kep')
    a = coord(1);
    e = coord(2);
    i = coord(3);
    OMEGA = coord(4);
    omega = coord(5);
    theta = coord(6);
    [r0, v0] = par2car(a, e, i, OMEGA, omega, theta, mu);
elseif strcmp(type, 'Car')
    r0 = coord(1:3);
    v0 = coord(4:6);
else 
    error('Not recognised. Use ''Kep'' or ''Car''.')
end

% orbit propagation and results
%[~, Y] = Orbit_Analysis(r0, v0, mu, t, 'non_perturbed', 'no_plot_orbit', 'no_plots');


% state vector
y0 = [r0; v0]; 
% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[~, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu), t, y0, options);    

r = Y(:, 1:3);
v = Y(:, 4:6);

x = r(:,1);
y = r(:,2);
z = r(:,3);
rnorm = vecnorm(r, 2, 2);

% calculation of alpha delta longitude and latitude
alpha = atan2(y,x);
delta = asin(z./rnorm);
thetaG = thetaG0 + omegaE*(t - t0);

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

end

