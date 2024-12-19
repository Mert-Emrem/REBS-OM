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
mu_Me       = astroConstants(11);                       %[km^3/s^2]     Mercury planetary costant
R_Me        = astroConstants(21);                       %[km]           Mercury mean radius
J2_Me       = astroConstants(31);                       %[-]            Mercury J2 gravitational harmonic coefficent
P_Me        = astroConstants(51);                       %[h]            Mercury sidereal rotation period
OMEGA_Me    = 2*pi/P_Me/60/60;                          %[rad/s]        Mercury angular velocity
eps_Me      = deg2rad(astroConstants(61));              %[rad]          Mercury axial tilt


% Orbit Keplerian elements
a_0         = 0.775e4;                                  %[km]       Initial semi-major axis
e_0         = 0.3528;                                   %[-]        Inital eccentricity
i_0         = deg2rad(18.52);                           %[rad]      Initial inclination
OM_0        = 0;                                        %           Can be changed by us
om_0        = pi/3;                                     %           Can be changed by us
theta_0     = 0;                                        %           Can be changed by us
k_el_0        = [a_0;e_0;i_0;OM_0;om_0;theta_0];        %           Keplerian elements' vector
T_orb       = 2*pi*sqrt(a_0^3/mu_Me);                   %[s]        Orbital period
start_date  = 0;                                        % date mjd2000

% Orbit initial conditions in cartesian coordinates
[rr_0, vv_0] = par2car (a_0, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_coord_0 = [rr_0; vv_0];

% Repeating GT ratio
m           = 3;                                        %[-]    rotations of the planet
k           = 5;                                        %[-]    revolutions of the satellite

% Options for the integrator
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% number of orbits and points for the simulations
n_orbits            = 1000;
n_points            = n_orbits * 250;
tspan = [0:100:n_orbits*T_orb];

%% Propagation and ground track of the unperturbed, assigned orbit

% Propagate orbit without perturbation acting on it
%[T_1,Y_1]           = orbit_propagator(rr_0 , vv_0 , mu_Me , 2 , options, n_points , n_orbits*T_orb);
[T_1,Y_1]            = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan, 'non_perturbed');

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

% modified semi-major axis
n = OMEGA_Me *k/m;
a_modif = (mu_Me/n^2)^(1/3);


% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
[rr_A_modif, vv_A_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);

% Orbit period with a_modif
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);
%%
% Propagate orbit with a_modif without perturbation acting on it
%[T_modif,Y_modif]                = orbit_propagator(rr_0_modif , vv_0_modif , mu_Me , 2 , options, n_points , n_orbits*T_orb_modif);
[T_modif,Y_modif]                 = Orbit_Analysis(rr_0_modif, vv_0_modif, mu_Me, tspan, 'non_perturbed');

% Calculate latitude and longitude of the ground track with a_modif
[alpha, delta, lon_modif, lat_modif]    = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);


% Convert latitude and longitude from radiants to degrees
lat_modif                       = lat_modif./pi .* 180;
lon_modif                       = lon_modif./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif, lat_modif, T_modif);
title('Ground Track of the unperturbated orbit (a_{modif})');
ylabel(colorbar,'Time');

%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator
% for the repeating ground track

tspan = [0:1000:1000*T_orb];

k_el_0_modif = [a_modif; k_el_0(2:6)];

tic
[T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me), mu_Me), tspan, k_el_0_modif, options);
elapsedtime_gauss_modif = toc;
fprintf ('Gaussian simulation for modified orbit: %.4f s. \n', elapsedtime_gauss_modif)

R_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);
V_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);

for j = 1:size(Y_pert_gauss_modif, 1)
    [rr, vv] = par2car(Y_pert_gauss_modif(j, 1), Y_pert_gauss_modif(j, 2), Y_pert_gauss_modif(j, 3), Y_pert_gauss_modif(j, 4), Y_pert_gauss_modif(j, 5), Y_pert_gauss_modif(j, 6), mu_Me);
    rr = rr';
    vv = vv';
    R_pert_gauss_modif(j,:) = rr;
    V_pert_gauss_modif(j,:) = vv;
end

% [R_pert_gauss, V_pert_gauss] = kep2car_vec(Y_pert_gauss(:, 1)', Y_pert_gauss(:, 2)', Y_pert_gauss(:, 3)', Y_pert_gauss(:, 4)', Y_pert_gauss(:, 5)', Y_pert_gauss(:, 6)', mu_Me)
% just with one orbit it saturates my RAM - the 3D Matrix dimension is too
% big

C_pert_gauss_modif = [R_pert_gauss_modif, V_pert_gauss_modif];

% Calculate latitude and longitude of the ground track with a_modif
[alpha, delta, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_pert_gauss_modif, OMEGA_Me, 0, 0, C_pert_gauss_modif);


% Convert latitude and longitude from radiants to degrees
lat_modif_pert                       = lat_modif_pert./pi .* 180;
lon_modif_pert                       = lon_modif_pert./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif_pert, lat_modif_pert, T_pert_gauss_modif);
title('Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
ylabel(colorbar,'Time');

%%
% Orbit Plot
Mercury_3D
hold on
scatter3(R_pert_gauss_modif(:,1), R_pert_gauss_modif(:,2), R_pert_gauss_modif(:,3), 3, tspan, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - gaussian propagator');

%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator

tspan = [0:100:n_orbits*T_orb];

tic
[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me), mu_Me), tspan, k_el_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation: %.4f s. \n', elapsedtime_gauss)

tic
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me), mu_Me), tspan, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation: %.4f s. \n', elapsedtime_cart)

%% Plot of perturbed orbit ground track from Gaussian propagator

R_pert_gauss = zeros(size(Y_pert_gauss,1), 3);
V_pert_gauss = zeros(size(Y_pert_gauss,1), 3);

for j = 1:size(Y_pert_gauss, 1)
    [rr, vv] = par2car(Y_pert_gauss(j, 1), Y_pert_gauss(j, 2), Y_pert_gauss(j, 3), Y_pert_gauss(j, 4), Y_pert_gauss(j, 5), Y_pert_gauss(j, 6), mu_Me);
    rr = rr';
    vv = vv';
    R_pert_gauss(j,:) = rr;
    V_pert_gauss(j,:) = vv;
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
Mercury_3D
hold on
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
Mercury_3D
hold on
%plot3(Y_pert_cart(:, 1), Y_pert_cart(:, 2), Y_pert_cart(:, 3));
scatter3(Y_pert_cart(:, 1), Y_pert_cart(:, 2), Y_pert_cart(:, 3), 3, tspan, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - cartesian propagator');

%% Error visualization between the two propagators

% keplerian parameters from the Cartesian propagator
for j = 1:size(Y_pert_cart, 1)
    [a_c, e_c, i_c, OM_c, om_c, th_c] = car2par(Y_pert_cart(j,1:3), Y_pert_cart(j,4:6), mu_Me);
    A(j) = a_c;
    E(j) = e_c;
    I(j) = i_c;
    RAAN(j) = OM_c;
    PER_AN(j) = om_c;
    TH(j) = th_c;
end

% the orbit propagators does not considered the angles in the [0, 2*pi]
% range authomatically while from car2par the angles are in the [0, 2*pi]
% range
for i = 1:size(Y_pert_gauss,1)
    while Y_pert_gauss(i,6)>2*pi
        Y_pert_gauss(i,6) = Y_pert_gauss(i,6) - 2*pi;
    end
end

% for i = 1:size(Y_pert_gauss,1)
%     while Y_pert_gauss(i,4)>2*pi
%         Y_pert_gauss(i,4) = Y_pert_gauss(i,4) - 2*pi;
%     end
% end
% 
% for i = 1:size(Y_pert_gauss,1)
%     while Y_pert_gauss(i,5)>2*pi
%         Y_pert_gauss(i,5) = Y_pert_gauss(i,5) - 2*pi;
%     end
% end

% semi-major axis 
figure
error_a = abs(A(:)-Y_pert_gauss(:,1))/a_0;
semilogy(tspan/T_orb, error_a, 'b');
grid on
title('a')

% eccentricity 
figure
error_e = abs(E(:)-Y_pert_gauss(:,2));
semilogy(tspan/T_orb, error_e, 'b');
grid on
title('e')

% inclination  
figure
error_i = abs(I(:)-Y_pert_gauss(:,3));
semilogy(tspan/T_orb, error_i, 'b');
grid on
title('i')

% RAAN  
figure
error_OM = abs(RAAN(:)-Y_pert_gauss(:,4))/(2*pi);
semilogy(tspan/T_orb, error_OM, 'b');
grid on
title('\Omega')

% pericenter anomaly
figure
error_om = abs(PER_AN(:)-Y_pert_gauss(:,5))/(2*pi);
semilogy(tspan/T_orb, error_om, 'b');
grid on
title('\omega')

% true anomaly
figure
error_th = abs(rad2deg(TH(:)-Y_pert_gauss(:,6)));
semilogy(tspan/T_orb, error_th, 'b');
grid on
title('\theta')

% error_x = (R_pert_gauss(:,1) - Y_pert_cart(:,1));
% error_y = (R_pert_gauss(:,2) - Y_pert_cart(:,2));
% error_z = (R_pert_gauss(:,3) - Y_pert_cart(:,3));
% 
% figure
% plot(tspan, error_x, 'r');
% hold on
% figure
% plot(tspan, error_y, 'b');
% figure
% plot(tspan, error_z, 'g');
% grid on
% %title('absolute value of the error on each component of the position vector')
% %legend('error_x', 'error_y', 'error_z');
% 
% error = vecnorm(R_pert_gauss, 2, 2) - vecnorm(Y_pert_cart(:, 1:3), 2,2);
% figure
% plot(tspan, error)
% grid on

%% Behaviour of keplerian parameters and filtering
figure
subplot(6,1,1)
a_filtered = movmean(Y_pert_gauss(:,1), T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,1));
hold on
plot(tspan./T_orb, a_filtered, 'k');
grid on                                                                                                   
title('a')


subplot(6,1,2)
e_filtered = movmean(Y_pert_gauss(:,2),  T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,2));
hold on
plot(tspan./T_orb, e_filtered, 'k');
grid on                                                                                                   
title('e')

subplot(6,1,3)
i_filtered = movmean(Y_pert_gauss(:,3),  T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,3));
hold on
plot(tspan./T_orb, i_filtered, 'k');
grid on                                                                                                   
title('i')

subplot(6,1,4)
OM_filtered = movmean(Y_pert_gauss(:,4),  T_orb);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,4)));
hold on
plot(tspan./T_orb, rad2deg(OM_filtered), 'k');
grid on                                                                                                   
title('\Omega')

subplot(6,1,5)
om_filtered = movmean(Y_pert_gauss(:,5),  T_orb);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,5)));
hold on
plot(tspan./T_orb, rad2deg(om_filtered), 'k');
grid on                                                                                                   
title('\omega')

subplot(6,1,6)
th_filtered = movmean(Y_pert_gauss(:,6),  T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,6));
hold on
plot(tspan./T_orb, th_filtered, 'k');
grid on                                                                                                   
title('\theta')

%%
figure
om_filtered = movmean(Y_pert_gauss(:,5),  T_orb);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,5)));
hold on
plot(tspan./T_orb, rad2deg(om_filtered), 'k');
grid on                                                                                                   
title('\omega')
%% 

figure

a_filtered = movmean(Y_pert_gauss(:,1), [T_orb/2 T_orb/2]);
plot(tspan./T_orb, Y_pert_gauss(:,1));
hold on
plot(tspan./T_orb, a_filtered, 'k');
grid on                                                                                                   
title('a')


figure
e_filtered = movmean(Y_pert_gauss(:,2),  [T_orb/2 T_orb/2]);
plot(tspan./T_orb, Y_pert_gauss(:,2));
hold on
plot(tspan./T_orb, e_filtered, 'k');
grid on                                                                                                   
title('e')

figure
i_filtered = movmean(Y_pert_gauss(:,3),  [T_orb/2 T_orb/2]);
plot(tspan./T_orb, Y_pert_gauss(:,3));
hold on
plot(tspan./T_orb, i_filtered, 'k');
grid on                                                                                                   
title('i')

figure
OM_filtered = movmean(Y_pert_gauss(:,4),  [T_orb/2 T_orb/2]);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,4)));
hold on
plot(tspan./T_orb, rad2deg(OM_filtered), 'k');
grid on                                                                                                   
title('\Omega')

figure
om_filtered = movmean(Y_pert_gauss(:,5),  [T_orb/2 T_orb/2]);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,5)));
hold on
plot(tspan./T_orb, rad2deg(om_filtered), 'k');
grid on                                                                                                   
title('\omega')

figure
th_filtered = movmean(Y_pert_gauss(:,6),  [T_orb/2 T_orb/2]);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,6)));
hold on
plot(tspan./T_orb, rad2deg(th_filtered), 'k');
grid on                                                                                                   
title('\theta')

%% 
ksun = astroConstants(4)
mu_Me
[ame, ~] = uplanet(0, 1)
rMe = 0.39 * astroConstants(2)
r_SOI = ame(1) * (mu_Me/ksun)^(2/5)
norm(rr_A_modif)





