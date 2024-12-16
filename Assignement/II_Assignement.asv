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
OM_0        = 0;                                %           Can be changed by us
om_0        = 0;                                %           Can be changed by us
theta_0     = 0;                                %           Can be changed by us
k_el        = [a_0;e_0;i_0;OM_0;om_0;theta_0];  %           Keplerian elements' vector
T_orb       = 2*pi*sqrt(a_0^3/mu_Me);           %[s]        Orbital period
start_date  = 0;

% Orbit initial conditions in cartesian coordinates
[rr_0, vv_0] = par2car (a_0, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

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
a_modif                         = a_groundTrack(k, m, mu_Me, OMEGA_Me);

% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);

% Propagate orbit with a_modif without perturbation acting on it
n_orbits                        = 1;
n_points                        = n_orbits * 500;
[T_modif,Y_modif]               = orbit_propagator(rr_0_modif , vv_0_modif , mu_Me , 2 , options, n_points , n_orbits*T_orb_modif);

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif, lat_modif]    = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);

% Convert latitude and longitude from radiants to degrees
lat_modif                       = lat_modif./pi .* 180;
lon_modif                       = lon_modif./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif, lat_modif, T_modif);
title('Ground Track of the unperturbated orbit (a_modif)');
ylabel(colorbar,'Time');

%% Introduce perturbations
[T_pert, Y_pert] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(T, k_el, J2_Me, mu_Me, R_E, day), mu_E, J2, R_Me, start_date), tspan, kep_0, options);






