%% Second assignement
%
% AIM:The PoliMi Space Agency wants to launch a Planetary Explorer Mission, 
% to perform planetary observations.
% As part of the mission analysis team, you are requested to carry out the 
% orbit analysis and ground track estimation. You have to study the effects
% of orbit perturbations, and compare different propagation methods. Also, 
% you have to characterize the ground track, and propose an orbit modifica-
% tion to reach a repeating ground track (for better observation of the 
% areas of interest).
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% N°GROUP: 2406 -REBS Project
%
% DATE: 
%
% -------------------------------------------------------------------------

clc
clear
close all

%% Inital data
% Central body: Mercury
% Perturbation: J2 + 3-body (Sun)

% Mercury physical properties
mu_Me       = astroConstants(11);               %[km^3/s^2]     Mercury planetary costant
R_Me        = astroConstants(21);               %[km]           Mercury mean radius
J2_Me       = astroConstants(31);               %[-]            Mercury J2 gravitational harmonic coefficent
P_Me        = astroConstants(51);               %[h]            Mercury sidereal rotation period
OMEGA_Me    = 2*pi/P_Me/60/60;                  %[rad/s]        Mercury angular velocity
eps_Me      = deg2rad(astroConstants(61));      %[rad]          Mercury axial tilt


% Orbit Keplerian elements
a_0         = 0.775e4;                          %[km]       Initial semi-major axis
e_0         = 0.3528;                           %[-]        Inital eccentricity
i_0         = deg2rad(18.52);                   %[rad]      Initial inclination
OM_0        = pi/8;                                %           Can be changed by us
om_0        = pi/3;                                %           Can be changed by us
theta_0     = 0;                                %           Can be changed by us
k_el_0        = [a_0;e_0;i_0;OM_0;om_0;theta_0];  %           Keplerian elements' vector
T_orb       = 2*pi*sqrt(a_0^3/mu_Me);           %[s]        Orbital period
start_date  = 0;

% Orbit initial conditions in cartesian coordinates
[rr_0, vv_0] = par2car (a_0, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_coord_0 = [rr_0; vv_0];

% Repeating GT ratio
m           = 3;                        %[-]    rotations of the planet
k           = 5;                        %[-]    revolutions of the satellite

% Options for the integrator
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%% Propagation and ground track of the unperturbed, assigned orbit

% Propagate orbit without perturbation acting on it
n_orbits            = 10;
n_points            = n_orbits * 500;
[T_1,Y_1]           = orbit_propagator(rr_0 , vv_0 , mu_Me , 2 , options, n_points , n_orbits*T_orb);

% Calculate latitude and longitude of the ground track
[~, ~, lon, lat]    = groundTrack_2(T_1, OMEGA_Me, 0, 0, Y_1);

% Convert latitude and longitude from radiants to degrees
lat                 = lat./pi .* 180;
lon                 = lon./pi .* 180;

% Plot ground track
plotGroundTrack(lon, lat, T_1);
title('Ground Track of the unperturbated, assigned orbit');
ylabel(colorbar,'Time');

%% Compute modified semi-major axis to obtain repeating orbits and groundtrak with a_modif
% NOTE: VERIFICARE VALIDITA'

% Calculate a_modif
%<<<<<<< Updated upstream
%a_modif                         = a_groundTrack(k, m, mu_Me, OMEGA_Me);
% =======
%a_modif                         = a_groundTrack(k, m, OMEGA_Me);
n = OMEGA_Me *k/m;
a_modif = (mu_Me/n^2)^(1/3);


% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
[rr_apo_modif, vv_apo_modif]    = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);

% Orbit period with a_modif
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);

% Propagate orbit with a_modif without perturbation acting on it
n_orbits                        = 10;
n_points                        = n_orbits * 500;
[T_modif,Y_modif]               = orbit_propagator(rr_0_modif , vv_0_modif , mu_Me , 2 , options, n_points , n_orbits*T_orb_modif);

% Calculate latitude and longitude of the ground track with a_modif
[alpha, delta, lon_modif, lat_modif]    = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);

figure
plot(lon_modif, lat_modif)

% Convert latitude and longitude from radiants to degrees
lat_modif                       = lat_modif./pi .* 180;
lon_modif                       = lon_modif./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif, lat_modif, T_modif);
title('Ground Track of the unperturbated orbit (a_{modif})');
ylabel(colorbar,'Time');

% =======
% figure
% plot(lon_modif, lat_modif)
% >>>>>>> Stashed changes
% 
% figure
% plot3(Y_modif(:,1), Y_modif(:,2), Y_modif(:,3));
% hold on
% plot3(0, 0, 0, LineWidth=70)
% grid on
% axis equal

%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator
clc
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
tspan = [0:1:10*T_orb];

[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date), mu_Me), tspan, k_el_0, options);
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date), mu_Me), tspan, cart_coord_0, options);

%% Plot of perturbed orbit ground track from Gaussian propagator

R_pert_gauss = zeros(size(Y_pert_gauss,1), 3);
V_pert_gauss = zeros(size(Y_pert_gauss,1), 3);

for i = 1:size(Y_pert_gauss, 1)
    [rr, vv] = par2car(Y_pert_gauss(i, 1), Y_pert_gauss(i, 2), Y_pert_gauss(i, 3), Y_pert_gauss(i, 4), Y_pert_gauss(i, 5), Y_pert_gauss(i, 6), mu_Me);
    rr = rr';
    vv = vv';
    R_pert_gauss(i,:) = rr;
    V_pert_gauss(i,:) = vv;
end

% [R_pert_gauss, V_pert_gauss] = kep2car_vec(Y_pert_gauss(:, 1)', Y_pert_gauss(:, 2)', Y_pert_gauss(:, 3)', Y_pert_gauss(:, 4)', Y_pert_gauss(:, 5)', Y_pert_gauss(:, 6)', mu_Me)
% just with one orbit it saturates my RAM - the 3D Matrix dimension is too
% big

C_pert_gauss = [R_pert_gauss, V_pert_gauss];
[alpha_pert, delta_pert, lon_pert_gauss, lat_pert_gauss]    = groundTrack_2(T_pert_gauss, OMEGA_Me, 0, 0, C_pert_gauss);

% Convert latitude and longitude from radiants to degrees
lat_pert_gauss                       = lat_pert_gauss./pi .* 180;
lon_pert_gauss                     = lon_pert_gauss./pi .* 180;

% Plot ground track
plotGroundTrack(lon_pert_gauss, lat_pert_gauss, T_pert_gauss);
title('Ground Track of the perturbated gauss orbit: J2 and 3BP-Sun');
ylabel(colorbar,'Time');

% Orbit Plot
figure
scatter3(R_pert_gauss(:,1), R_pert_gauss(:,2), R_pert_gauss(:,3), 3, tspan, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - gaussian propagator');

%% Plot of perturbed orbit ground track from Cartesian propagator

[~, ~, lon_pert_cart, lat_pert_cart]    = groundTrack_2(T_pert_cart, OMEGA_Me, 0, 0, Y_pert_cart);

% Convert latitude and longitude from radiants to degrees
lat_pert_cart                       = lat_pert_cart./pi .* 180;
lon_pert_cart                       = lon_pert_cart./pi .* 180;

% Plot ground track
plotGroundTrack(lon_pert_cart, lat_pert_cart, T_pert_cart);
title('Ground Track of the perturbated gauss orbit: J2 and 3BP-Sun');
ylabel(colorbar,'Time');

% Orbit Plot 
figure
scatter3(Y_pert_cart(:, 1), Y_pert_cart(:, 2), Y_pert_cart(:, 3), 3, tspan, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - cartesian propagator');

%% Error visualization between the two propagators

error_x = abs(R_pert_gauss(:,1) - Y_pert_cart(:,1));
error_y = abs(R_pert_gauss(:,2) - Y_pert_cart(:,2));
error_z = abs(R_pert_gauss(:,3) - Y_pert_cart(:,3));

figure
plot(tspan, error_x, 'r');
hold on
plot(tspan, error_y, 'b');
plot(tspan, error_z, 'g');
grid on
title('absolute value of the error on each component of the position vector')
legend('error_x', 'error_y', 'error_z');
error = vecnorm(R_pert_gauss, 2, 2) - vecnorm(Y_pert_cart(:, 1:3), 2,2);
figure
plot(tspan, error)
grid on
