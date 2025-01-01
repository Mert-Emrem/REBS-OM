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

set(0,'defaultfigurecolor',[1 1 1])

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
r_mean_Me   = 5.791e+7; %0.39 * astroConstants(2)                 %[km]           Mercury's orbit mean ratio
mu_Sun      = astroConstants(4);                        %[km^3/s^2]     Sun planetary costant  

% Orbit Keplerian elements
a_0         = 0.775e4;                                  %[km]       Initial semi-major axis
e_0         = 0.3528;                                   %[-]        Inital eccentricity
i_0         = deg2rad(18.52);                           %[rad]      Initial inclination
OM_0        = deg2rad(150);                                     %           Can be changed by us
om_0        = pi/3;                                     %           Can be changed by us
theta_0     = 0;                                        %           Can be changed by us
k_el_0      = [a_0;e_0;i_0;OM_0;om_0;theta_0];          %           Keplerian elements' vector
T_orb       = 2*pi*sqrt(a_0^3/mu_Me);                   %[s]        Orbital period
start_date  = date2mjd2000([2041 3 28 0 0 0]);          %[mjd2000]  Starting date of the planetary mission                           % date mjd2000

% Orbit initial conditions in cartesian coordinates
[rr_0, vv_0] = par2car (a_0, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_coord_0 = [rr_0; vv_0];

% Repeating GT ratio
m           = 3;                                        %[-]    rotations of the planet
k           = 5;                                        %[-]    revolutions of the satellite

% Options for the integrator
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% number of orbits and points for the simulations
n_orbits            = 5;
n_points            = n_orbits * 700;
tspan = [0:10:n_orbits*T_orb];

%% Propagation and ground track of the unperturbed assigned orbit

% Propagate orbit without perturbation acting on it
[T_1,Y_1]            = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan, 'non_perturbed');
delta_lambda = rad2deg(T_orb * OMEGA_Me);
fprintf('GT: Delta Lambda = %f \n', delta_lambda);

% Calculate latitude and longitude of the ground track
[~, ~, lon, lat]    = groundTrack_2(T_1, OMEGA_Me, 0, 0, Y_1);

% Plot ground track
plotGroundTrack(lon, lat, T_1);
title('Ground Track of the unperturbated assigned orbit');
ylabel(colorbar,'Time');


N = [cos(OM_0); sin(OM_0); 0];
h = cross(rr_0, vv_0)/norm(cross(rr_0, vv_0));
y = cross(h, N);
e = rotation_vector_Rodrigues(N, h, om_0);
p = cross(h, e);

Mercury_3D
hold on
plot3(Y_1(:,1), Y_1(:,2), Y_1(:,3), 'b', LineWidth=2);
quiver3(0, 0, 0, 1.5*R_Me, 0, 0, 'Color',[178/255, 0, 0], LineWidth=1); 
text(1.6*R_Me, 0, 0, 'X', 'FontSize', 9, 'Color', [178/255, 0, 0]);
quiver3(0, 0, 0, 0, 1.5*R_Me, 0, 'Color', [178/255, 0, 0], LineWidth=1); 
text(0, 1.6*R_Me,  0, 'Y', 'FontSize', 9, 'Color', [178/255, 0, 0]);
quiver3(0, 0, 0, 0, 0, 1.5*R_Me, 'Color',[178/255, 0, 0], LineWidth=1); 
text(0, 0, 1.6*R_Me,  'Z', 'FontSize', 9, 'Color', [178/255, 0, 0]);
quiver3(0, 0, 0, 1.5*R_Me*N(1), 1.5*R_Me*N(2), N(3), 'Color', [255/255, 165/255, 0/255], LineWidth=1); 
text(1.5*R_Me*cos(OM_0), 1.5*R_Me*sin(OM_0), 0, 'N', 'FontSize', 9, 'Color', [255/255, 165/255, 0/255]);
quiver3(0, 0, 0, 1.5*R_Me*h(1), 1.5*R_Me*h(2), 1.5*R_Me*h(3), 'k', LineWidth=1); 
text( 1.6*R_Me*h(1), 1.6*R_Me*h(2), 1.6*R_Me*h(3), 'h', 'FontSize', 9, 'Color', 'k');
quiver3(0, 0, 0, 1.5*R_Me*e(1), 1.5*R_Me*e(2), 1.5*R_Me*e(3), 'k', LineWidth=1); 
text(1.6*R_Me*e(1), 1.6*R_Me*e(2), 1.6*R_Me*e(3), 'e', 'FontSize', 9, 'Color', 'k');
quiver3(0, 0, 0, 1.5*R_Me*p(1), 1.5*R_Me*p(2), 1.5*R_Me*p(3), 'k', LineWidth=1); 
text(1.6*R_Me*p(1), 1.6*R_Me*p(2), 1.6*R_Me*p(3), 'p', 'FontSize', 9, 'Color', 'k');
grid on
axis equal
title('Nominal Orbit and reference frames');
legend('', 'Nominal Orbit', 'Equatorial Reference Frame', '', '', 'Ascending Node', 'Orbital Reference Frame', '', '');

%% Compute modified semi-major axis to obtain repeating orbits and groundtrack with a_modif

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

% Plot ground track
plotGroundTrack(lon_modif, lat_modif, T_modif);
title('Ground Track of the unperturbated orbit (a_{modif})');
ylabel(colorbar,'Time');

% PLot Orbit
Mercury_3D
hold on
plot3(Y_1(:,1), Y_1(:,2), Y_1(:,3), 'c', LineWidth=2);
plot3(Y_modif(:,1), Y_modif(:,2), Y_modif(:,3), 'b', LineWidth=2);
grid on
axis equal
title('Unperturbed orbit for the repeating ground track');
legend('Mercury', 'Nominal Orbit', 'Modified orbit for RGT')

%% Check: is the apocenter of the modified orbit inside the sphere of interest of Mercury?

[rr_A_modif, vv_A_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);
 
r_SOI = r_mean_Me * (mu_Me/mu_Sun)^(2/5);
r_apo = norm(rr_A_modif);
n_perc = 0.5;
limit = n_perc*r_SOI/(1+e_0);

a_rep = zeros(30,30);
RATIO = zeros(30,30);

flag = 1;
if r_apo > r_SOI 
    fprintf('Apocenter radius is larger than the sphere of interest of Mercury\n');
    for K = [1:1:30]
        for M = [1:1:30]
            k_m_ratio = K/M;
            a_modif_new = (mu_Me*(M/(K*OMEGA_Me))^2)^(1/3);
            a_rep(K,M) = a_modif_new;
            RATIO(K,M) = k_m_ratio;
            if a_modif_new <= ((n_perc*r_SOI)/(1+e_0))
                a_rep_valid(flag) = a_modif_new;
                RATIO_valid(flag) = k_m_ratio;
                k_valid(flag) = K;
                m_valid(flag) = M;
                flag = flag+1;   
            end
        end
    end
else 
    fprintf('Apocenter radius is smaller than the sphere of interest of Mercury\n');
end

a_modif_new = a_rep_valid(1);
[r_apo_new, ~] = par2car(a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
r_apo_new_norm = norm(r_apo_new);
perc = r_apo_new_norm/r_SOI;
fprintf('For the new repeating GT, the condition of r_apo <= %.2f *r_SOI is respected: r_apo_new/r_SOI = %f \n', n_perc, perc);
k_modif_new = k_valid(1);
m_modif_new = m_valid(1);
fprintf('Valid value of semi-major axix: %f \n', a_modif_new);
fprintf('Valid value of k/m: %.1f, where k = %.1f and m = %.1f \n', k_modif_new/m_modif_new, k_modif_new, m_modif_new);


% Orbit initial conditions with a_modif_new in cartesian coordinates
[rr_0_modif_new, vv_0_modif_new]        = par2car (a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif_new
T_orb_modif_new                 = 2*pi*sqrt(a_modif_new^3/mu_Me);
fprintf('Valid value of period: %.2f \n', T_orb_modif_new);

tspan_2 = [0:100:k_modif_new*T_orb_modif_new];
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
plot3(Y_1(:,1), Y_1(:,2), Y_1(:,3), 'Color', [34/255, 139/255, 34/255], LineWidth=2);
plot3(Y_modif(:,1), Y_modif(:,2), Y_modif(:,3), 'b', LineWidth=2);
plot3(rr_A_modif(1), rr_A_modif(2), rr_A_modif(3), 'k*', LineWidth=1)
plot3(Y_modif_new(:,1), Y_modif_new(:,2), Y_modif_new(:,3), 'r', LineWidth=2);
grid on
axis equal
title('Unperturbed orbit for the repeating ground track');
[X_sphere,Y_sphere,Z_sphere] = sphere(50);
X_sphere = X_sphere*r_SOI;
Y_sphere = Y_sphere*r_SOI;
Z_sphere = Z_sphere*r_SOI;
SOI = surf(X_sphere,Y_sphere,Z_sphere);
set(SOI, 'FaceColor', 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'c', 'EdgeAlpha', 0.2); 
legend('Mercury', 'Nominal Orbit','Repeating orbit', 'Apocenter', 'New repeating orbit', 'SOI');

%% TO BE DELETED

% if r_apo > r_SOI 
%     fprintf('Apocenter radius is larger than the sphere of interest of Mercury\n');
%     n_perc = 0.9;
%     r_apo_new = n_perc*r_SOI;
%     a_modif_new = r_apo_new / (1+e_0);
%     k_m_ratio_new = OMEGA_Me * sqrt((a_modif_new^3)/mu_Me);
%     k_m_ratio_new = 1/round(k_m_ratio_new, 4);
%     k_new = 5;
%     m_new = k_new/k_m_ratio_new;
%     m_new = floor(m_new);
% %     m_new = 3;
% %     k_new = m_new * k_m_ratio_new;
% %     k_new = floor(k_new);
%     a_modif_new = (mu_Me*(m_new/(k_new*OMEGA_Me))^2)^(1/3);
%     r_apo_new = a_modif_new * (1+e_0);
% else 
%     fprintf('Apocenter radius is smaller than the sphere of interest of Mercury\n');
% end
% %
% Orbit initial conditions with a_modif_new in cartesian coordinates
% [rr_0_modif_new, vv_0_modif_new]        = par2car (a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
% 
% Orbit period with a_modif_new
% T_orb_modif_new                 = 2*pi*sqrt(a_modif_new^3/mu_Me);
% 
% tspan_2 = [0:100:k_modif_new*T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
% [T_modif_new,Y_modif_new]                 = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_2, 'non_perturbed');
% 
% Calculate latitude and longitude of the ground track with a_modif
% [~, ~, lon_modif_new, lat_modif_new]    = groundTrack_2(T_modif_new, OMEGA_Me, 0, 0, Y_modif_new);
% 
% Plot ground track
% plotGroundTrack(lon_modif_new, lat_modif_new, T_modif_new);
% title('Ground Track of the unperturbated orbit (a_{modif,new})');
% ylabel(colorbar,'Time');
% 
% Orbit plot
% Mercury_3D
% hold on
% plot3(Y_modif(:,1), Y_modif(:,2), Y_modif(:,3), 'b', LineWidth=2);
% plot3(rr_A_modif(1), rr_A_modif(2), rr_A_modif(3), 'r*', LineWidth=1)
% plot3(Y_modif_new(:,1), Y_modif_new(:,2), Y_modif_new(:,3), 'm', LineWidth=2);
% grid on
% axis equal
% title('Unperturbed orbit for the repeating ground track');
% [X_sphere,Y_sphere,Z_sphere] = sphere(50);
% X_sphere = X_sphere*r_SOI;
% Y_sphere = Y_sphere*r_SOI;
% Z_sphere = Z_sphere*r_SOI;
% SOI = surf(X_sphere,Y_sphere,Z_sphere);
% set(SOI, 'FaceColor', 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'c', 'EdgeAlpha', 0.2); 
% legend('Mercury', 'Repeating orbit', 'Apocenter', 'New repeating orbit', 'SOI');

%% Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
T_orb_modif                    = 2*pi*sqrt(a_modif^3/mu_Me);
tspan_1 = [0:100:T_orb_modif];
% SIDE NOTES: T_orb_modif è 105 volte più grande dell'orbita assegnata,
% quindi integrare per molte orbite crea errori di convergenza causa bassa
% toleranza

k_el_0_modif = [a_modif; e_0; i_0; OM_0; om_0; theta_0];

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

% Plot ground track
plotGroundTrack(lon_modif_pert, lat_modif_pert, T_pert_gauss_modif);
title('Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
ylabel(colorbar,'Time');

% %Orbit Plot
% Mercury_3D
% hold on
% scatter3(R_pert_gauss_modif(:,1), R_pert_gauss_modif(:,2), R_pert_gauss_modif(:,3), 3, tspan_1, 'filled');
% colormap('jet')
% colorbar
% grid on
% axis equal
% title('orbit - gaussian propagator');


%% CORRECTED Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
tspan_1 = [0:100:10*T_orb_modif_new];

k_el_0_modif = [a_modif_new; e_0; i_0; OM_0;om_0;theta_0];
[rmod0, vmod0] = par2car(a_modif_new, e_0, i_0, OM_0,om_0,theta_0, mu_Me);

cart_0_mod = [rmod0; vmod0];
tic
[T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, k_el_0_modif, options);
elapsedtime_gauss_modif = toc;
fprintf ('Gaussian simulation for modified orbit: %.4f s. \n', elapsedtime_gauss_modif)

tic
[T_pert_gauss_modifc, Y_pert_gauss_modifc] = ode113(@(t, cc) eq_motion_cartesian(t, cc, @(t, cc) acc_pert_function_cartesian(t, cc, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, cart_0_mod, options);
elapsedtime_cart_modif = toc;
fprintf ('Gaussian simulation for modified orbit: %.4f s. \n', elapsedtime_cart_modif)

R_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);
V_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);

for j = 1:size(Y_pert_gauss_modif, 1)
    [rr, vv] = par2car(Y_pert_gauss_modif(j, 1), Y_pert_gauss_modif(j, 2), Y_pert_gauss_modif(j, 3), Y_pert_gauss_modif(j, 4), Y_pert_gauss_modif(j, 5), Y_pert_gauss_modif(j, 6), mu_Me);
    rr = rr';
    vv = vv';
    R_pert_gauss_modif(j,:) = rr;
    V_pert_gauss_modif(j,:) = vv;
end

C_pert_gauss_modif = [R_pert_gauss_modif, V_pert_gauss_modif];

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_pert_gauss_modif, OMEGA_Me, 0, 0, C_pert_gauss_modif);


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

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pertc, lat_modif_pertc]    = groundTrack_2(T_pert_gauss_modifc, OMEGA_Me, 0, 0, Y_pert_gauss_modifc);


% Plot ground track
plotGroundTrack(lon_modif_pertc, lat_modif_pertc, T_pert_gauss_modifc);
title('Ground Track of the perturbated orbit cart (a_{modif}) - J2 and 3BP');
ylabel(colorbar,'Time');

%Orbit Plot
Mercury_3D
hold on
scatter3(Y_pert_gauss_modifc(:,1), Y_pert_gauss_modifc(:,2), Y_pert_gauss_modifc(:,3), 3, tspan_1, 'filled');
colormap('jet')
colorbar
grid on
axis equal
title('orbit - cartesian propagator');


%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator

tspan = [0:50:500*T_orb];

tic
[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, k_el_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation: %.4f s. \n', elapsedtime_gauss)

tic
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation: %.4f s. \n', elapsedtime_cart);

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
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';
grid on
axis equal
title('orbit - gaussian propagator');

%% Plot of perturbed orbit ground track from Cartesian propagator

[~, ~, lon_pert_cart, lat_pert_cart]    = groundTrack_2(T_pert_cart, OMEGA_Me, 0, 0, Y_pert_cart);

% Plot ground track
plotGroundTrack(lon_pert_cart, lat_pert_cart, T_pert_cart);
title('Ground Track of the perturbated gauss orbit: J2 and 3BP-Sun');
ylabel(colorbar,'Time');

% Orbit Plot 
Mercury_3D
hold on
scatter3(Y_pert_cart(:, 1), Y_pert_cart(:, 2), Y_pert_cart(:, 3), 3, tspan, 'filled');
colormap
colorbar
c = colorbar;
c.Label.String = 'Time [s]';
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

%% 7TH REQUEST
clear all
close all
clc

% clear all
% close all
% clc
% 
% data_open = fopen('DecDec.txt', 'r');
% 
% 
% Leggi i dati come stringhe di caratteri
% textData = textscan(data_open, '%s', 'Delimiter', ';');
% 
% Chiudi il file
% fclose(data_open);
% 
% remove spaces
% cleanedData = strrep(textData{1}, ' ', '');
% 
% ephemerides = [];
%
% % strings to double
% for i = 1:length(cleanedData)
%     numArray = str2double(strsplit(cleanedData{i}, ','));
%     ephemerides = [ephemerides; numArray];
% end

% OBJECTS CHOSEN:
% - LEO ORBIT: Envisat
% - MEO ORBIT: -------
%%
% Data for Earth as primary planet
mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 12 14 0 0 0]);

%% ENVISAT
% upload the information from envisat_eph
load('envisateph.mat');

a_eph_envi = envisateph(:,10);
e_eph_envi = envisateph(:,1);
i_eph_envi = deg2rad(envisateph(:,3));
OM_eph_envi = deg2rad(envisateph(:,4));
om_eph_envi = deg2rad(envisateph(:,5));
th_eph_envi = deg2rad(envisateph(:,9));

T_orb_eph_envi = 2*pi*sqrt(a_eph_envi(1)^3/mu_E);
nn = 864000/T_orb_eph_envi;
tspan3 = [0:60:10*24*60*60];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

eph_0_envi = [a_eph_envi(1); e_eph_envi(1); i_eph_envi(1); OM_eph_envi(1); om_eph_envi(1); th_eph_envi(1)];

for i = 1:length(a_eph_envi)
    [rr, vv] = par2car(a_eph_envi(i), e_eph_envi(i), i_eph_envi(i), OM_eph_envi(i), om_eph_envi(i), th_eph_envi(i), mu_E);
    r(i, :) = rr;
    v(i, :) = vv;
end
Terra_3D
hold on
plot3(r(:,1), r(:,2), r(:,3))
grid on 
axis equal

tic
[T_pert_eph_envi, Y_pert_eph_envi] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_envi, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph_envi, 1)
    [rr, vv] = par2car(Y_pert_eph_envi(i,1), Y_pert_eph_envi(i,2), Y_pert_eph_envi(i,3), Y_pert_eph_envi(i,4), Y_pert_eph_envi(i,5), Y_pert_eph_envi(i,6), mu_E);
    rrr(i, :) = rr;
    vvv(i, :) = vv;
end

hold on
scatter3(rrr(:,1), rrr(:,2), rrr(:,3), 3, T_pert_eph_envi, 'filled');
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';

%% Behaviour

for i = 1:size(Y_pert_eph_envi,1)
    while Y_pert_eph_envi(i,6)>=2*pi
        Y_pert_eph_envi(i,6) = Y_pert_eph_envi(i,6) - 2*pi;
    end
end

for i = 1:size(Y_pert_eph_envi,1)
    while Y_pert_eph_envi(i,4)>=2*pi
        Y_pert_eph_envi(i,4) = Y_pert_eph_envi(i,4) - 2*pi;
    end
     while Y_pert_eph_envi(i,4)<0
        Y_pert_eph_envi(i,4) = Y_pert_eph_envi(i,4) + 2*pi;
    end
end

for i = 1:size(Y_pert_eph_envi,1)
    while Y_pert_eph_envi(i,5)>=2*pi
        Y_pert_eph_envi(i,5) = Y_pert_eph_envi(i,5) - 2*pi;
    end
    while Y_pert_eph_envi(i,5)<0
        Y_pert_eph_envi(i,5) = Y_pert_eph_envi(i,5) + 2*pi;
    end
end

figure
vect = tspan3;
plot(vect, a_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,1));
grid on
title('a')
legend('ephemerides', 'our propagator')

figure
plot(vect, e_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,2));
grid on
title('e')
legend('ephemerides', 'our propagator')

figure
plot(vect, i_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,3));
grid on
title('i')
legend('ephemerides', 'our propagator')

figure
plot(vect, OM_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,4));
grid on
title('\Omega')
legend('ephemerides', 'our propagator')

figure
plot(vect, om_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,5));
grid on
title('\omega')
legend('ephemerides', 'our propagator')

figure
plot(vect, th_eph_envi(:));
hold on
plot(vect, Y_pert_eph_envi(:,6));
grid on
title('\theta')
legend('ephemerides', 'our propagator')


%% Errors - not required for the assignment
% for i = 1:size(Y_pert_eph,1)
%     while Y_pert_eph(i,6)>=2*pi
%         Y_pert_eph(i,6) = Y_pert_eph(i,6) - 2*pi;
%     end
% end
% 
% for i = 1:size(Y_pert_eph,1)
%     while Y_pert_eph(i,4)>=2*pi
%         Y_pert_eph(i,4) = Y_pert_eph(i,4) - 2*pi;
%     end
%      while Y_pert_eph(i,4)<0
%         Y_pert_eph(i,4) = Y_pert_eph(i,4) + 2*pi;
%     end
% end
% 
% for i = 1:size(Y_pert_eph,1)
%     while Y_pert_eph(i,5)>=2*pi
%         Y_pert_eph(i,5) = Y_pert_eph(i,5) - 2*pi;
%     end
%     while Y_pert_eph(i,5)<0
%         Y_pert_eph(i,5) = Y_pert_eph(i,5) + 2*pi;
%     end
% end
% 
% % semi-major axis 
% figure
% error_a = abs(a_eph(:)-Y_pert_eph(:,1))/a_eph(1);
% semilogy(tspan3/T_orb_eph, error_a, 'b');
% grid on
% title('a')
% 
% % eccentricity 
% figure
% error_e = abs(e_eph(:)-Y_pert_eph(:,2));
% semilogy(tspan3/T_orb_eph, error_e, 'b');
% grid on
% title('e')
% 
% % inclination  
% figure
% error_i = abs(i_eph(:)-Y_pert_eph(:,3));
% semilogy(tspan3/T_orb_eph, error_i, 'b');
% grid on
% title('i')
% 
% % RAAN  
% figure
% error_OM = abs(OM_eph(:)-Y_pert_eph(:,4))/(2*pi);
% semilogy(tspan3/T_orb_eph, error_OM, 'b');
% grid on
% title('\Omega')
% 
% % pericenter anomaly
% figure
% error_om = abs(om_eph(:)-Y_pert_eph(:,5))/(2*pi);
% semilogy(tspan3/T_orb_eph, error_om, 'b');
% grid on
% title('\omega')
% 
% % true anomaly
% figure
% error_th = abs(rad2deg(th_eph(:)-Y_pert_eph(:,6))); %./abs(Y_pert_gauss(:,6));
% semilogy(tspan3/T_orb_eph, error_th, 'b');
% grid on
% title('\theta')


%% GOES 13
% upload the information from goes13_eph
load('goes13eph.mat');

a_eph_goes = goes13eph(:,10);
e_eph_goes = goes13eph(:,1);
i_eph_goes = deg2rad(goes13eph(:,3));
OM_eph_goes = deg2rad(goes13eph(:,4));
om_eph_goes = deg2rad(goes13eph(:,5));
th_eph_goes = deg2rad(goes13eph(:,9));

T_orb_eph_goes = 2*pi*sqrt(a_eph_goes(1)^3/mu_E);
nn = 864000/T_orb_eph_goes;
tspan3 = [0:60:10*24*60*60];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

eph_0_goes = [a_eph_goes(1); e_eph_goes(1); i_eph_goes(1); OM_eph_goes(1); om_eph_goes(1); th_eph_goes(1)];

for i = 1:length(a_eph_goes)
    [rr, vv] = par2car(a_eph_goes(i), e_eph_goes(i), i_eph_goes(i), OM_eph_goes(i), om_eph_goes(i), th_eph_goes(i), mu_E);
    r_goes(i, :) = rr;
    v_goes(i, :) = vv;
end

Terra_3D
hold on
plot3(r_goes(:,1), r_goes(:,2), r_goes(:,3))
grid on 
axis equal

tic
[T_pert_eph_goes, Y_pert_eph_goes] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_goes, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph_goes, 1)
    [rr, vv] = par2car(Y_pert_eph_goes(i,1), Y_pert_eph_goes(i,2), Y_pert_eph_goes(i,3), Y_pert_eph_goes(i,4), Y_pert_eph_goes(i,5), Y_pert_eph_goes(i,6), mu_E);
    rrr_goes(i, :) = rr;
    vvv_goes(i, :) = vv;
end

hold on
scatter3(rrr_goes(:,1), rrr_goes(:,2), rrr_goes(:,3), 3, T_pert_eph_goes, 'filled');
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';

%% Behaviour

% the orbit propagators does not considered the angles in the [0, 2*pi]
% range authomatically while from car2par the angles are in the [0, 2*pi]
% range
for i = 1:size(Y_pert_eph_goes,1)
    while Y_pert_eph_goes(i,6)>=2*pi
        Y_pert_eph_goes(i,6) = Y_pert_eph_goes(i,6) - 2*pi;
    end
end

for i = 1:size(Y_pert_eph_goes,1)
    while Y_pert_eph_goes(i,4)>=2*pi
        Y_pert_eph_goes(i,4) = Y_pert_eph_goes(i,4) - 2*pi;
    end
     while Y_pert_eph_goes(i,4)<0
        Y_pert_eph_goes(i,4) = Y_pert_eph_goes(i,4) + 2*pi;
    end
end

for i = 1:size(Y_pert_eph_goes,1)
    while Y_pert_eph_goes(i,5)>=2*pi
        Y_pert_eph_goes(i,5) = Y_pert_eph_goes(i,5) - 2*pi;
    end
    while Y_pert_eph_goes(i,5)<0
        Y_pert_eph_goes(i,5) = Y_pert_eph_goes(i,5) + 2*pi;
    end
end


figure
vect = tspan3;
plot(vect, a_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,1));
grid on
title('a')
legend('ephemerides', 'our propagator')

figure
plot(vect, e_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,2));
grid on
title('e')
legend('ephemerides', 'our propagator')

figure
plot(vect, i_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,3));
grid on
title('i')
legend('ephemerides', 'our propagator')

figure
plot(vect, OM_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,4));
grid on
title('\Omega')
legend('ephemerides', 'our propagator')

figure
plot(vect, om_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,5));
grid on
title('\omega')
legend('ephemerides', 'our propagator')

figure
plot(vect, th_eph_goes(:));
hold on
plot(vect, Y_pert_eph_goes(:,6));
grid on
title('\theta')
legend('ephemerides', 'our propagator')


