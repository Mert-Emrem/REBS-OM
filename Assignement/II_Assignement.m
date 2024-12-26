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

addpath 'Functions'\
addpath 'Functions'\timeConversion\time\


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
r_mean_Me   = 0.39 * astroConstants(2);                 %[km]           Mercury's orbit mean ratio
mu_Sun      = astroConstants(4);                        %[km^3/s^2]     Sun planetary costant  

% Orbit Keplerian elements
a_0         = 0.775e4;                                  %[km]       Initial semi-major axis
e_0         = 0.3528;                                   %[-]        Inital eccentricity
i_0         = deg2rad(18.52);                           %[rad]      Initial inclination
OM_0        = pi/6;                                     %           Can be changed by us
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
n_orbits            = 100;
n_points            = n_orbits * 500;
tspan = [0:10:n_orbits*T_orb];

%% Propagation and ground track of the unperturbed assigned orbit

% Propagate orbit without perturbation acting on it
%[T_1,Y_1]           = orbit_propagator(rr_0 , vv_0 , mu_Me , 2 , options, n_points , n_orbits*T_orb);
[T_1,Y_1]            = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan, 'non_perturbed');

% Calculate latitude and longitude of the ground track
[~, ~, lon, lat]    = groundTrack_2(T_1, OMEGA_Me, 0, 0, Y_1);

% % Convert latitude and longitude from radiants to degrees
% lat                 = lat./pi .* 180;
% lon                 = lon./pi .* 180;

% Plot ground track
plotGroundTrack(lon, lat, T_1);
title('Ground Track of the unperturbated assigned orbit');
ylabel(colorbar,'Time');

%% Compute modified semi-major axis to obtain repeating orbits and groundtrak with a_modif

% modified semi-major axis
n = OMEGA_Me *k/m;
a_modif = (mu_Me/n^2)^(1/3);

% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);

% after k S/C orbits the ground track will repeat equal to itself
tspan_1 = [0:100:k*T_orb_modif];

% Propagate orbit with a_modif without perturbation acting on it
[T_modif,Y_modif]                 = Orbit_Analysis(rr_0_modif, vv_0_modif, mu_Me, tspan_1, 'non_perturbed');

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif, lat_modif]    = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);

% % Convert latitude and longitude from radiants to degrees
% lat_modif                       = lat_modif./pi .* 180;
% lon_modif                       = lon_modif./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif, lat_modif, T_modif);
title('Ground Track of the unperturbated orbit (a_{modif})');
ylabel(colorbar,'Time');


%% Check: is the apocenter of the modified orbit inside the sphere of interest of Mercury?

[rr_A_modif, vv_A_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);
 
r_SOI = r_mean_Me * (mu_Me/mu_Sun)^(2/5);
r_apo = norm(rr_A_modif);
if r_apo > r_SOI 
    fprintf('Apocenter radius is larger than the sphere of interest of Mercury\n');
    n_perc = 0.9;
    r_apo_new = n_perc*r_SOI;
    a_modif_new = r_apo_new / (1+e_0);
    k_m_ratio_new = OMEGA_Me * sqrt((a_modif_new^3)/mu_Me);
    %k_m_ratio_new = round(k_m_ratio_new, 4);
    k_new = 5;
    m_new = k_new/k_m_ratio_new;
    m_new = floor(m_new);
%     m_new = 3;
%     k_new = m_new * k_m_ratio_new;
%     k_new = floor(k_new);
%     a_modif_new = (mu_Me*(k_new/(m_new*OMEGA_Me))^2)^(1/3)
%     r_apo_new = a_modif_new * (1+e_0);
else 
    fprintf('Apocenter radius is smaller than the sphere of interest of Mercury\n');
end

% Orbit initial conditions with a_modif_new in cartesian coordinates
[rr_0_modif_new, vv_0_modif_new]        = par2car (a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif_new
T_orb_modif_new                 = 2*pi*sqrt(a_modif_new^3/mu_Me);

tspan_2 = [0:100:713*T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
[T_modif_new,Y_modif_new]                 = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_2, 'non_perturbed');

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_new, lat_modif_new]    = groundTrack_2(T_modif_new, OMEGA_Me, 0, 0, Y_modif_new);

% Plot ground track
plotGroundTrack(lon_modif_new, lat_modif_new, T_modif_new);
title('Ground Track of the unperturbated orbit (a_{modif,new})');
ylabel(colorbar,'Time');

% Orbit plot
Mercury_3D
hold on
plot3(Y_modif(:,1), Y_modif(:,2), Y_modif(:,3), 'b', LineWidth=2);
plot3(rr_A_modif(1), rr_A_modif(2), rr_A_modif(3), 'r*', LineWidth=1)
plot3(Y_modif_new(:,1), Y_modif_new(:,2), Y_modif_new(:,3), 'm', LineWidth=2);
grid on
axis equal
title('Unperturbed orbit for the repeating ground track');
[X_sphere,Y_sphere,Z_sphere] = sphere(50);
X_sphere = X_sphere*r_SOI;
Y_sphere = Y_sphere*r_SOI;
Z_sphere = Z_sphere*r_SOI;
SOI = surf(X_sphere,Y_sphere,Z_sphere);
set(SOI, 'FaceColor', 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'c', 'EdgeAlpha', 0.2); 
legend('Mercury', 'Repeating orbit', 'Apocenter', 'New repeating orbit', 'SOI');

%% Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);
tspan_1 = [0:1000:T_orb_modif];
% SIDE NOTES: T_orb_modif è 105 volte più grande dell'orbita assegnata,
% quindi integrare per molte orbite crea errori di convergenza causa bassa
% toleranza

k_el_0_modif = [a_modif; k_el_0(2:6)];

tic
[T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, k_el_0_modif, options);
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
[~, ~, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_pert_gauss_modif, OMEGA_Me, 0, 0, C_pert_gauss_modif);

%  % Convert latitude and longitude from radiants to degrees
% lat_modif_pert                       = lat_modif_pert./pi .* 180;
% lon_modif_pert                       = lon_modif_pert./pi .* 180;

% Plot ground track
plotGroundTrack(lon_modif_pert, lat_modif_pert, T_pert_gauss_modif);
title('Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
ylabel(colorbar,'Time');

%Orbit Plot
Mercury_3D
hold on
scatter3(R_pert_gauss_modif(:,1), R_pert_gauss_modif(:,2), R_pert_gauss_modif(:,3), 3, tspan_1, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - gaussian propagator');

%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator

%tspan = [0:100:n_orbits*T_orb];

tic
[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, k_el_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation: %.4f s. \n', elapsedtime_gauss)

tic
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, cart_coord_0, options);
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

% % Convert latitude and longitude from radiants to degrees
% lat_pert_cart                       = lat_pert_cart./pi .* 180;
% lon_pert_cart                       = lon_pert_cart./pi .* 180;

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
    while Y_pert_gauss(i,6)>=2*pi
        Y_pert_gauss(i,6) = Y_pert_gauss(i,6) - 2*pi;
    end
end

for i = 1:size(Y_pert_gauss,1)
    while Y_pert_gauss(i,4)>=2*pi
        Y_pert_gauss(i,4) = Y_pert_gauss(i,4) - 2*pi;
    end
     while Y_pert_gauss(i,4)<0
        Y_pert_gauss(i,4) = Y_pert_gauss(i,4) + 2*pi;
    end
end

for i = 1:size(Y_pert_gauss,1)
    while Y_pert_gauss(i,5)>=2*pi
        Y_pert_gauss(i,5) = Y_pert_gauss(i,5) - 2*pi;
    end
    while Y_pert_gauss(i,5)<0
        Y_pert_gauss(i,5) = Y_pert_gauss(i,5) + 2*pi;
    end
end

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
error_th = abs(rad2deg(TH(:)-Y_pert_gauss(:,6))); %./abs(Y_pert_gauss(:,6));
semilogy(tspan/T_orb, error_th, 'b');
grid on
title('\theta')

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

%% SIDE NOTES
% - with only J2 effect active we see, as espected from theory:
%      - regression of OMEGA
%      - procession of omega 
% - with J2+3BP effects active we see:
%      - regression of OMEGA
%      - regression of omega 
%      - procession of i

%%
clear all
close all
clc

data_open = fopen('DecDec.txt', 'r');

% Leggi i dati come stringhe di caratteri
textData = textscan(data_open, '%s', 'Delimiter', ';');

% Chiudi il file
fclose(data_open);

% remove spaces
cleanedData = strrep(textData{1}, ' ', '');

ephemerides = [];

% strings to double
for i = 1:length(cleanedData)
    numArray = str2double(strsplit(cleanedData{i}, ','));
    ephemerides = [ephemerides; numArray];
end

mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 12 14 0 0 0]);


a_eph = ephemerides(:,10);
e_eph = ephemerides(:,1);
i_eph = deg2rad(ephemerides(:,3));
OM_eph = deg2rad(ephemerides(:,4));
om_eph = deg2rad(ephemerides(:,5));
th_eph = deg2rad(ephemerides(:,9));

T_orb_eph = 2*pi*sqrt(a_eph(1)^3/mu_E);
nn = 864000/T_orb_eph;
tspan3 = [0:6000:nn*T_orb_eph];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

eph_0 = [a_eph(1); e_eph(1); i_eph(1); OM_eph(1); om_eph(1); th_eph(1)];

for i = 1:145
    [rr, vv] = par2car(a_eph(i), e_eph(i), i_eph(i), OM_eph(i), om_eph(i), th_eph(i), mu_E);
    r(i, :) = rr;
    v(i, :) = vv;
end

plot3(r(:,1), r(:,2), r(:,3))
grid on 
axis equal

tic
[T_pert_eph, Y_pert_eph] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph, 1)
    [rr, vv] = par2car(Y_pert_eph(i,1), Y_pert_eph(i,2), Y_pert_eph(i,3), Y_pert_eph(i,4), Y_pert_eph(i,5), Y_pert_eph(i,6), mu_E);
    rrr(i, :) = rr;
    vvv(i, :) = vv;
end

hold on
scatter3(rrr(:,1), rrr(:,2), rrr(:,3), 3, T_pert_eph, 'filled');
colormap('jet')
colorbar
%% Behaviour

figure
vect = linspace(0, T_orb_eph, 145);
plot(vect, a_eph(:));
hold on
plot(vect, Y_pert_eph(:,1));
grid on
title('a')


figure
plot(vect, e_eph(:));
hold on
plot(vect, Y_pert_eph(:,2));
grid on
title('e')

figure
plot(vect, i_eph(:));
hold on
plot(vect, Y_pert_eph(:,3));
grid on
title('i')

figure
plot(vect, OM_eph(:));
hold on
plot(vect, Y_pert_eph(:,4));
grid on
title('\Omega')

figure
plot(vect, om_eph(:));
hold on
plot(vect, Y_pert_eph(:,5));
grid on
title('\omega')

figure
plot(vect, th_eph(:));
hold on
plot(vect, Y_pert_eph(:,6));
grid on
title('\theta')


%% Errors
for i = 1:size(Y_pert_eph,1)
    while Y_pert_eph(i,6)>=2*pi
        Y_pert_eph(i,6) = Y_pert_eph(i,6) - 2*pi;
    end
end

for i = 1:size(Y_pert_eph,1)
    while Y_pert_eph(i,4)>=2*pi
        Y_pert_eph(i,4) = Y_pert_eph(i,4) - 2*pi;
    end
     while Y_pert_eph(i,4)<0
        Y_pert_eph(i,4) = Y_pert_eph(i,4) + 2*pi;
    end
end

for i = 1:size(Y_pert_eph,1)
    while Y_pert_eph(i,5)>=2*pi
        Y_pert_eph(i,5) = Y_pert_eph(i,5) - 2*pi;
    end
    while Y_pert_eph(i,5)<0
        Y_pert_eph(i,5) = Y_pert_eph(i,5) + 2*pi;
    end
end

% semi-major axis 
figure
error_a = abs(a_eph(:)-Y_pert_eph(:,1))/a_eph(1);
semilogy(tspan3/T_orb_eph, error_a, 'b');
grid on
title('a')

% eccentricity 
figure
error_e = abs(e_eph(:)-Y_pert_eph(:,2));
semilogy(tspan3/T_orb_eph, error_e, 'b');
grid on
title('e')

% inclination  
figure
error_i = abs(i_eph(:)-Y_pert_eph(:,3));
semilogy(tspan3/T_orb_eph, error_i, 'b');
grid on
title('i')

% RAAN  
figure
error_OM = abs(OM_eph(:)-Y_pert_eph(:,4))/(2*pi);
semilogy(tspan3/T_orb_eph, error_OM, 'b');
grid on
title('\Omega')

% pericenter anomaly
figure
error_om = abs(om_eph(:)-Y_pert_eph(:,5))/(2*pi);
semilogy(tspan3/T_orb_eph, error_om, 'b');
grid on
title('\omega')

% true anomaly
figure
error_th = abs(rad2deg(th_eph(:)-Y_pert_eph(:,6))); %./abs(Y_pert_gauss(:,6));
semilogy(tspan3/T_orb_eph, error_th, 'b');
grid on
title('\theta')




