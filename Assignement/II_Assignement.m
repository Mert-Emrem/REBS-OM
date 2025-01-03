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
a_0         = 0.775e4;                                  %[km]       Initial semi-major axis (assigned)
e_0         = 0.3528;                                   %[-]        Inital eccentricity (assigned)
i_0         = deg2rad(18.52);                           %[rad]      Initial inclination (assigned)
OM_0        = deg2rad(150);                             %[rad]      Initial RAAN (chosen)
om_0        = pi/3;                                     %[rad]      Initial pericenter anomaly (chosen)
theta_0     = 0;                                        %[rad]      Initial true anomaly (chosen)
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
n_orbits            = 500;
%n_points            = n_orbits * 700;
span                = T_orb/500;
tspan = [0:span:n_orbits*T_orb];


%% Propagation and ground track of the unperturbed assigned orbit

% 1 orbit
tspan = [0:span/10:1*T_orb];
% Propagate orbit without perturbation acting on it
[T_1,Y_1] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan, 'non_perturbed');
delta_lambda = rad2deg(T_orb * OMEGA_Me);
fprintf('GT: Delta Lambda = %f \n', delta_lambda);

% Calculate latitude and longitude of the ground track
[~, ~, lon, lat] = groundTrack_2(T_1, OMEGA_Me, 0, 0, Y_1);

% 1 days
tspan_1_day = [0:span/10:24*60*60];
% Propagate orbit without perturbation acting on it
[T_10,Y_10] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan_1_day, 'non_perturbed');
% Calculate latitude and longitude of the ground track
[~, ~, lon_10, lat_10] = groundTrack_2(T_10, OMEGA_Me, 0, 0, Y_10);

% 10 days
tspan_10_days = [0:span/10:10*24*60*60];
% Propagate orbit without perturbation acting on it
[T_100,Y_100] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan_10_days, 'non_perturbed');
% Calculate latitude and longitude of the ground track
[~, ~, lon_100, lat_100] = groundTrack_2(T_100, OMEGA_Me, 0, 0, Y_100);

% Plot ground track
figure
subplot(3,1,1)
plotGroundTrack(gca, lon, lat, T_1);
title('Unperturbated Ground Track: 1 period');
subplot(3,1,2)
plotGroundTrack(gca, lon_10, lat_10, T_10);
title('Unperturbated Ground Track: 1 day');
subplot(3,1,3)
plotGroundTrack(gca, lon_100, lat_100, T_100);
title('Unperturbated Ground Track: 10 days');

% reference system vectors for plot
N = [cos(OM_0); sin(OM_0); 0];
h = cross(rr_0, vv_0)/norm(cross(rr_0, vv_0));
y = cross(h, N);
e = rotation_vector_Rodrigues(N, h, om_0);
p = cross(h, e);

% Nomonal orbit plot with reference systems
Mercury_3D
hold on
plotOrbit(Y_1, 2);
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
legend('Mercury', 'Nominal Orbit da fun' ,'Equatorial Reference Frame', '', '', 'Ascending Node', 'Orbital Reference Frame', '', '');


%% Compute modified semi-major axis to obtain repeating orbits and groundtrack with a_modif

% modified semi-major axis
n = OMEGA_Me *k/m;
a_modif = (mu_Me/n^2)^(1/3);

% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif]        = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif
T_orb_modif                     = 2*pi*sqrt(a_modif^3/mu_Me);
delta_lambda_modif = rad2deg(T_orb_modif * OMEGA_Me);
fprintf('Repeating GT: Delta Lambda  = %f \n', delta_lambda_modif);

% after k S/C orbits the ground track will repeat equal to itself
tspan_1 = [0:span:k*T_orb_modif];

% Propagate orbit with a_modif without perturbation acting on it
[T_modif,Y_modif] = Orbit_Analysis(rr_0_modif, vv_0_modif, mu_Me, tspan_1, 'non_perturbed');

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif, lat_modif] = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);

% Plot ground track
figure
plotGroundTrack(gca, lon_modif, lat_modif, T_modif);
title('Unperturbed Repeating Ground Track, k=5, m=3');

% Plot Orbit
Mercury_3D
hold on
plotOrbit(Y_1, 2);
plotOrbit(Y_modif, 2);
grid on
axis equal
title('Unperturbed orbit for the repeating ground track');
legend('Mercury', 'Nominal Orbit', 'Modified orbit for RGT')


%% Check: is the apocenter of the modified orbit inside the sphere of interest of Mercury?

% apocenter and pericenter position vector
[rr_A_modif, ~]        = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);
r_apo = norm(rr_A_modif);
[rr_P_modif, ~]        = par2car (a_modif, e_0, i_0, OM_0, om_0, 0, mu_Me);
r_peri = norm(rr_P_modif);

r_SOI = r_mean_Me * (mu_Me/mu_Sun)^(2/5);        %sphere of influence of Mercury
n_perc = 0.5;                                    % limit on the percentage of SOI for apocenter

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
k_new = k_valid(1);
m_modif_new = m_valid(1);
fprintf('Valid value of semi-major axix: %f \n', a_modif_new);
fprintf('Valid value of k/m: %.1f, where k = %.1f and m = %.1f \n', k_new/m_modif_new, k_new, m_modif_new);

% Orbit initial conditions with a_modif_new in cartesian coordinates
[rr_0_modif_new, vv_0_modif_new] = par2car (a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif_new
T_orb_modif_new  = 2*pi*sqrt(a_modif_new^3/mu_Me);
fprintf('Valid value of period: %.2f \n', T_orb_modif_new);
delta_lambda_new = rad2deg(T_orb_modif_new * OMEGA_Me);
fprintf('New Repating GT: Delta Lambda = %.2f \n\n', delta_lambda_new);

tspan_2 = [0:span:T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
[T_modif_new,Y_modif_new] = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_2, 'non_perturbed');
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_new, lat_modif_new] = groundTrack_2(T_modif_new, OMEGA_Me, 0, 0, Y_modif_new);

tspan_15 = [0:span:k_new*T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
[T_modif_new_15,Y_modif_new_15] = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_15, 'non_perturbed');
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_new_15, lat_modif_new_15] = groundTrack_2(T_modif_new_15, OMEGA_Me, 0, 0, Y_modif_new_15);

% Plot ground track
figure
subplot(2,1,1)
plotGroundTrack(gca, lon_modif_new, lat_modif_new, T_modif_new);
title('Unperturbed Repeating Ground Track, k=15, m=1: 1 orbit');
subplot(2,1,2)
plotGroundTrack(gca, lon_modif_new_15, lat_modif_new_15, T_modif_new_15);
title('Unperturbed Repeating Ground Track, k=15, m=1: 15 orbits');

% Orbit plot: nominal, repeating and new repeating and SOI
Mercury_3D
hold on
plotOrbit(Y_1, 2);
plotOrbit(Y_modif, 2);
plot3(rr_A_modif(1), rr_A_modif(2), rr_A_modif(3), 'k*', LineWidth=1)
plotOrbit(Y_modif_new, 2);
plotSOI(r_SOI)
grid on
axis equal
title('Unperturbed orbit for the repeating ground track and the new one');
legend('Mercury', 'Nominal Orbit','Repeating orbit', 'Apocenter', 'New repeating orbit', 'SOI');

% Orbit plot: nominal and new repeating
Mercury_3D
hold on
plotOrbit(Y_1, 2)
plotOrbit(Y_modif_new, 2)
grid on
axis equal
title('Unperturbed orbit for the new repeating ground track'); 
legend('Mercury', 'Nominal Orbit', 'New repeating orbit');


%%
% Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator, Cartesian, Orbit_Analysis
% %% Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator
% 
% options = odeset('RelTol',1e-13,'AbsTol',1e-14);                                 % SIDE NOTES: T_orb_modif è 105 volte più grande dell'orbita assegnata,
% T_orb_modif                    = 2*pi*sqrt(a_modif^3/mu_Me);                   % quindi integrare per molte orbite crea errori di convergenza causa bassa
% time_check = T_orb_modif/T_orb;                                                % toleranza
% tspan_1 = [0:span:50*T_orb_modif];                                             % it encounter some sort of discontinuity 
% 
% % to approximately see where the sun is
% [mer_pos, k_sun] = uplanet(start_date, 1);
% mer_pos = norm(mer_pos);
% [rr_mer, vv_sun] = par2car(mer_pos, 0, -eps_Me, 0, 0, 0, k_sun);
% rr_sun = -rr_mer; 
% 
% % propagation with Gauss
% k_el_0_modif = [a_modif; e_0; i_0; OM_0; om_0; theta_0];
% [rr_0_test, vv_0_test] = par2car(a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
% 
% tic
% [T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, k_el_0_modif, options);
% elapsedtime_gauss_modif = toc;
% fprintf ('Gaussian simulation for modified orbit: %.4f s. \n', elapsedtime_gauss_modif)
% 
% % to see value of the perturbation
% vect_gt_mod_big = acc_pert_function_test_value(T_pert_gauss_modif,  Y_pert_gauss_modif, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1);
% 
% figure
% subplot(6, 2, 1)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,1));
% title('a')
% subplot(6, 2, 2)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,2));
% title('e')
% subplot(6, 2, 3)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,3));
% title('i')
% subplot(6, 2, 4)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,4));
% title('OM')
% subplot(6, 2, 5)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,5));
% title('om')
% subplot(6, 2, 6)
% plot(T_pert_gauss_modif, Y_pert_gauss_modif(:,6));
% title('theta')
% 
% [T_test,Y_test]            = Orbit_Analysis(rr_0_test, vv_0_test, mu_Me, tspan, 'perturbed', J2_Me, R_Me);
% % Calculate latitude and longitude of the ground track with a_modif
% [~, ~, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_test, OMEGA_Me, 0, 0, Y_test);
% 
% % Plot ground track
% plotGroundTrack(lon_modif_pert, lat_modif_pert, T_test);
% title('ORBIT ANALYSIS Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
% ylabel(colorbar,'Time');
% 
% R_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);
% V_pert_gauss_modif = zeros(size(Y_pert_gauss_modif,1), 3);
% 
% for j = 1:size(Y_pert_gauss_modif, 1)
%     [rr, vv] = par2car(Y_pert_gauss_modif(j, 1), Y_pert_gauss_modif(j, 2), Y_pert_gauss_modif(j, 3), Y_pert_gauss_modif(j, 4), Y_pert_gauss_modif(j, 5), Y_pert_gauss_modif(j, 6), mu_Me);
%     rr = rr';
%     vv = vv';
%     R_pert_gauss_modif(j,:) = rr;
%     V_pert_gauss_modif(j,:) = vv;
% end
% 
% C_pert_gauss_modif = [R_pert_gauss_modif, V_pert_gauss_modif];
% 
% % Calculate latitude and longitude of the ground track with a_modif
% [~, ~, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_pert_gauss_modif, OMEGA_Me, 0, 0, C_pert_gauss_modif);
% 
% % Plot ground track
% plotGroundTrack(lon_modif_pert, lat_modif_pert, T_pert_gauss_modif);
% title('GAUSS Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
% ylabel(colorbar,'Time');
% 
% % Orbit Plot 
% Mercury_3D
% hold on
% plot3(R_pert_gauss_modif(:,1), R_pert_gauss_modif(:,2), R_pert_gauss_modif(:,3))
% %plot3(rr_sun(1), rr_sun(2), rr_sun(3), 'y*', LineWidth=15);
% grid on
% axis equal
% title('orbit - gaussian propagator');
% 
% % cartesian
% 
% % cartesia initial condition
% [rr_0_test, vv_0_test] = par2car(a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
% c_0_modif = [rr_0_test; vv_0_test];
% 
% % cartesian propagator
% tic
% [T_pert_cart_modif, Y_pert_cart_modif] = ode113(@(t, k_el) eq_motion_cartesian(t, k_el, @(t, k_el) acc_pert_function_cartesian(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, c_0_modif, options);
% elapsedtime_cart_modif = toc;
% fprintf ('Cartesian simulation for modified orbit: %.4f s. \n', elapsedtime_cart_modif)
% 
% vect_a_c = acc_pert_function_cartesian_test_value(T_pert_cart_modif, Y_pert_cart_modif, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1);
% vect_a_c_norm = vecnorm(vect_a_c, 2, 2);
% acc_grav = - mu_Me./(vecnorm(Y_pert_cart_modif(:, 1:3), 2, 2)).^3 .* Y_pert_cart_modif(:, 1:3);
% acc_grav_norm = vecnorm(acc_grav, 2, 2);
% 
% figure
% title('Pertrubed nominal orbit');
% subplot(2,1,1)
% plot(T_pert_cart_modif, vect_a_c_norm)
% title('perturbing acceleration')
% subplot(2,1,2)
% plot(T_pert_cart_modif, acc_grav_norm)
% title('gravitational acceleration')
% grid on
% 
% 
% % tic
% % [T_pert_cart_modif, Y_pert_cart_modif] = ode113(@(t, k_el) eq_motion_cartesian(t, k_el, @(t, k_el) acc_pert_function_cartesian(t, k_el, astroConstants(33), astroConstants(13), astroConstants(23), start_date, astroConstants(63), 3), mu_Me), tspan_1, c_0_modif, options);
% % elapsedtime_cart_modif = toc;
% % fprintf ('Cartesian simulation for modified orbit: %.4f s. \n', elapsedtime_cart_modif)
% 
% figure
% subplot(6, 2, 1)
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,1));
% title('a')       
% subplot(6, 2, 2)           
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,2));
% title('e')
% subplot(6, 2, 3)
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,3));
% title('i')
% subplot(6, 2, 4)
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,4));
% title('OM')
% subplot(6, 2, 5)
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,5));
% title('om')
% subplot(6, 2, 6)
% plot(T_pert_cart_modif, Y_pert_cart_modif(:,6));
% title('theta')
% 
% % % Calculate latitude and longitude of the ground track with a_modif
% % [~, ~, lon_modif_pert, lat_modif_pert]    = groundTrack_2(T_pert_cart_modif, OMEGA_Me, 0, 0, Y_pert_cart_modif);
% % 
% % % Plot ground track
% % plotGroundTrack(lon_modif_pert, lat_modif_pert, T_pert_cart_modif);
% % title('CARTESIAN Ground Track of the perturbated orbit cartesian (a_{modif}) - J2 and 3BP');
% % ylabel(colorbar,'Time');
% 
% 
% % Orbit Plot 
% Mercury_3D
% hold on
% plot3(Y_pert_cart_modif(:,1), Y_pert_cart_modif(:,2), Y_pert_cart_modif(:,3))
% scatterOrbit(Y_pert_cart_modif, tspan_1)
% %plot3(rr_sun(1), rr_sun(2), rr_sun(3), 'y*', LineWidth=15);
% r_SOI_earth = 149598000 * (astroConstants(13)/k_sun)^(2/5);
% %plotSOI(r_SOI_earth)
% grid on
% axis equal
% title('orbit - cartesian propagator');

%% CORRECTED Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator

% new time span
tspan_1 = [0:span:k_new*T_orb_modif_new];

% set of initial conditions with the new semi-major axis
k_el_0_modif = [a_modif_new; e_0; i_0; OM_0; om_0; theta_0];
[r_mod_0, v_mod_0] = par2car(a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_0_modif = [r_mod_0; v_mod_0];

% propagation with Gauss Planetary equations
tic
[T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, k_el_0_modif, options);
elapsedtime_gauss_modif = toc;
fprintf ('Gaussian simulation for modified repeating GT orbit: %.4f s. \n', elapsedtime_gauss_modif)

% check value of acceleration
v_a_mod = acc_pert_function_test_value(T_pert_gauss_modif, Y_pert_gauss_modif, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1);

% from Keplerian parameters to cartesian coordinates
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
figure
plotGroundTrack(gca, lon_modif_pert, lat_modif_pert, T_pert_gauss_modif);
title('Ground Track of the perturbated orbit (a_{modif}) - J2 and 3BP');
ylabel(colorbar,'Time');

%Orbit Plot
Mercury_3D
hold on
scatterOrbit(R_pert_gauss_modif, tspan_1);
grid on
axis equal
title('orbit - gaussian propagator');

% Propagation with Cartesian Equation of Motion 
% Computational time: highere than cartesian, but it does not require the
% conversion from keplerian to cartesian for the ground track plot

% 1 orbit
tspan_1 = [0:span:T_orb_modif_new];
tic
[T_pert_cart_modif, Y_pert_cart_modif] = ode113(@(t, cc) eq_motion_cartesian(t, cc, @(t, cc) acc_pert_function_cartesian(t, cc, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, cart_0_modif, options);
elapsedtime_cart_modif = toc;
fprintf ('Cartesian simulation for modified repeating GT orbit: %.4f s. \n\n', elapsedtime_cart_modif);
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pertc, lat_modif_pertc]    = groundTrack_2(T_pert_cart_modif, OMEGA_Me, 0, 0, Y_pert_cart_modif);

% 15 orbits
tspan_1 = [0:span:k_new*T_orb_modif_new];
tic
[T_pert_cart_modif_15, Y_pert_cart_modif_15] = ode113(@(t, cc) eq_motion_cartesian(t, cc, @(t, cc) acc_pert_function_cartesian(t, cc, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, cart_0_modif, options);
elapsedtime_cart_modif = toc;
fprintf ('Cartesian simulation for modified repeating GT orbit: %.4f s. \n\n', elapsedtime_cart_modif);
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pertc_15, lat_modif_pertc_15]    = groundTrack_2(T_pert_cart_modif_15, OMEGA_Me, 0, 0, Y_pert_cart_modif_15);

% Plot ground track
figure
subplot(2,2,1)
plotGroundTrack(gca, lon_modif_new, lat_modif_new, T_modif_new);
title('Unperturbed Repeating Ground Track, k=15, m=1: 1 orbit');
subplot(2,2,3)
plotGroundTrack(gca, lon_modif_new_15, lat_modif_new_15, T_modif_new_15);
title('Unperturbed Repeating Ground Track, k=15, m=1: 15 orbits');
subplot(2,2,2)
plotGroundTrack(gca, lon_modif_pertc, lat_modif_pertc, T_pert_cart_modif);
title('Perturbed Repeating Ground Track, k=15, m=1: 1 orbit');
subplot(2,2,4)
plotGroundTrack(gca, lon_modif_pertc_15, lat_modif_pertc_15, T_pert_cart_modif_15);
title('Perturbed Repeating Ground Track, k=15, m=1: 15 orbits');

% Orbit Plot
Mercury_3D
hold on
scatterOrbit(Y_pert_cart_modif, tspan_1);
grid on
axis equal
title('orbit - cartesian propagator');


%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator

% Nominal Initial Conditions
k_el_0      = [a_0;e_0;i_0;OM_0;om_0;theta_0];  

% new tspan appropriate to see the perturbations for the 3D plot
tspan = [0:span:1000*T_orb];

% Gaussian Propagator
tic
[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, k_el_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation for 100 orbits: %.4f s. \n', elapsedtime_gauss)

%%
% % check valute of acceleration in cart ref
% vect_a_g = acc_pert_function_test_value(T_pert_gauss, Y_pert_gauss, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1);
% vect_a_g_norm = vecnorm(vect_a_g, 2, 2);
% 
% figure
% plot(T_pert_gauss, vect_a_g_norm)
% grid on

% retrieve the matrix of the position vector and of the velocity vector
R_pert_gauss = zeros(size(Y_pert_gauss,1), 3);
for j = 1:size(Y_pert_gauss, 1)
    [rr, ~] = par2car(Y_pert_gauss(j, 1), Y_pert_gauss(j, 2), Y_pert_gauss(j, 3), Y_pert_gauss(j, 4), Y_pert_gauss(j, 5), Y_pert_gauss(j, 6), mu_Me);
    rr = rr';
    R_pert_gauss(j,:) = rr;
end

% Plot Orbit
Mercury_3D
hold on
scatterOrbit(R_pert_gauss, tspan);
grid on
axis equal
title('Perturbed orbit over 100 periods');

%%
% Cartesian Propagator
tic
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation over 100 orbits: %.4f s. \n', elapsedtime_cart);

% check valute of acceleration in cart ref
% vect_a_c = acc_pert_function_cartesian_test_value(T_pert_cart, Y_pert_cart, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1);
% vect_a_c_norm = vecnorm(vect_a_c, 2, 2);
% acc_grav = - mu_Me./(vecnorm(Y_pert_cart(:, 1:3), 2, 2)).^3 .* Y_pert_cart(:, 1:3);
% acc_grav_norm = vecnorm(acc_grav, 2, 2);

% figure
% title('Pertrubed nominal orbit');
% subplot(2,1,1)
% plot(T_pert_cart, vect_a_c_norm)
% title('perturbing acceleration')
% subplot(2,1,2)
% plot(T_pert_cart, acc_grav_norm)
% title('gravitational acceleration')
% grid on
%%
% Orbit Plot Cartesian
Mercury_3D
hold on
scatterOrbit(Y_pert_cart, tspan);
grid on
axis equal
title('orbit revolution over 100 periods - cartesian propagator');


%% Plot of perturbed orbit ground track from Gaussian and Cartesian propagators

% appropriate number of orbits to analize keplerian parameters behaviour
n_orb_GT = 100;

% more appropriate tspan for the GT
tspan_GT = [0:span/10:n_orb_GT*T_orb];

tic
[T_pert_gauss_GT, Y_pert_gauss_GT] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT, k_el_0, options);
elapsedtime_gauss_GT = toc;
fprintf ('Gaussian simulation: %.4f s. \n\n', elapsedtime_gauss_GT)


%%
% retrieve the matrix of the position vector and of the velocity vector
R_pert_gauss_GT = zeros(size(Y_pert_gauss_GT,1), 3);
V_pert_gauss_GT = zeros(size(Y_pert_gauss_GT,1), 3);
for j = 1:size(Y_pert_gauss_GT, 1)
    [rr, vv] = par2car(Y_pert_gauss_GT(j, 1), Y_pert_gauss_GT(j, 2), Y_pert_gauss_GT(j, 3), Y_pert_gauss_GT(j, 4), Y_pert_gauss_GT(j, 5), Y_pert_gauss_GT(j, 6), mu_Me);
    rr = rr';
    vv = vv';
    R_pert_gauss_GT(j,:) = rr;
    V_pert_gauss_GT(j,:) = vv;
end
C_pert_gauss_GT = [R_pert_gauss_GT, V_pert_gauss_GT];
[~, ~, lon_pert_gauss_GT, lat_pert_gauss_GT] = groundTrack_2(T_pert_gauss_GT, OMEGA_Me, 0, 0, C_pert_gauss_GT);

% Plot ground track
figure
plotGroundTrack(gca, lon_pert_gauss_GT, lat_pert_gauss_GT, T_pert_gauss_GT);
title('Ground Track of the perturbated orbit (Gauss) for 5 orbits: J2 and 3BP-Sun');


%% Cartesian Propagator

% 1 orbit
tspan_GT = [0:span/10:T_orb];
tic
[T_pert_cart_GT, Y_pert_cart_GT] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT, lat_pert_cart_GT] = groundTrack_2(T_pert_cart_GT, OMEGA_Me, 0, 0, Y_pert_cart_GT);

% 1 day
tspan_GT_1d = [0:span/10:24*60*60];
tic
[T_pert_cart_GT_1d, Y_pert_cart_GT_1d] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT_1d, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT_1d, lat_pert_cart_GT_1d] = groundTrack_2(T_pert_cart_GT_1d, OMEGA_Me, 0, 0, Y_pert_cart_GT_1d);

% 10 days
tspan_GT_10d = [0:span/10:10*24*60*60];
tic
[T_pert_cart_GT_10d, Y_pert_cart_GT_10d] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT_10d, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT_10d, lat_pert_cart_GT_10d] = groundTrack_2(T_pert_cart_GT_10d, OMEGA_Me, 0, 0, Y_pert_cart_GT_10d);


%%

% Plot ground track
figure
subplot(3,2,1)
plotGroundTrack(gca, lon, lat, T_1);
title('Unperturbated Ground Track: 1 period');
subplot(3,2,3)
plotGroundTrack(gca, lon_10, lat_10, T_10);
title('Unperturbated Ground Track: 1 day');
subplot(3,2,5)
plotGroundTrack(gca, lon_100, lat_100, T_100);
title('Unperturbated Ground Track: 10 days');
subplot(3,2,2)
plotGroundTrack(gca,lon_pert_cart_GT, lat_pert_cart_GT, T_pert_cart_GT);
title('Ground Track of the perturbated orbit (Cart) for 1 orbit');
subplot(3,2,4)
plotGroundTrack(gca,lon_pert_cart_GT_1d, lat_pert_cart_GT_1d, T_pert_cart_GT_1d);
title('Ground Track of the perturbated orbit (Cart) for 1 day');
subplot(3,2,6)
plotGroundTrack(gca,lon_pert_cart_GT_10d, lat_pert_cart_GT_10d, T_pert_cart_GT_10d);
title('Ground Track of the perturbated orbit (Cart) for 10 days');


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

% unwrap of the vectors
RAAN = unwrap(RAAN);
PER_AN = unwrap(PER_AN);
TH = unwrap(TH);
%%
%Y_pert_gauss(:,6) = unwrap(Y_pert_gauss(:,6));
%
% the orbit propagators does not considered the angles in the [0, 2*pi]
% range authomatically while from car2par the angles are in the [0, 2*pi]
% range
% for i = 1:size(Y_pert_gauss,1)
%     while Y_pert_gauss(i,6)>=2*pi
%         Y_pert_gauss(i,6) = Y_pert_gauss(i,6) - 2*pi;
%     end
% end
% 
% for i = 1:size(Y_pert_gauss,1)
%     while Y_pert_gauss(i,4)>=2*pi
%         Y_pert_gauss(i,4) = Y_pert_gauss(i,4) - 2*pi;
%     end
%      while Y_pert_gauss(i,4)<0
%         Y_pert_gauss(i,4) = Y_pert_gauss(i,4) + 2*pi;
%     end
% end
% 
% for i = 1:size(Y_pert_gauss,1)
%     while Y_pert_gauss(i,5)>=2*pi
%         Y_pert_gauss(i,5) = Y_pert_gauss(i,5) - 2*pi;
%     end
%     while Y_pert_gauss(i,5)<0
%         Y_pert_gauss(i,5) = Y_pert_gauss(i,5) + 2*pi;
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
error_th = abs(rad2deg(TH(:)-Y_pert_gauss(:,6)))./abs(rad2deg(Y_pert_gauss(:,6)));
semilogy(tspan/T_orb, error_th, 'b');
grid on
title('\theta')

%% Subplot figure

% semi-major axis 
figure
subplot(3,2,1)
error_a = abs(A(:)-Y_pert_gauss(:,1))/a_0;
semilogy(tspan/T_orb, error_a, 'b');
grid on
title('a')
ylabel('|a_{cart} - a_{gauss}|/a_0')
xlabel('Time in numbers of periods')

% eccentricity 
subplot(3,2,2)
error_e = abs(E(:)-Y_pert_gauss(:,2));
semilogy(tspan/T_orb, error_e, 'b');
grid on
title('e')
ylabel('|e_{cart} - e_{gauss}|')
xlabel('Time in numbers of periods')

% inclination  
subplot(3,2,3)
error_i = abs(I(:)-Y_pert_gauss(:,3));
semilogy(tspan/T_orb, error_i, 'b');
grid on
title('i')
ylabel('|i_{cart} - i_{gauss}|')
xlabel('Time in numbers of periods')

% RAAN  
subplot(3,2,4)
error_OM = abs(RAAN(:)-Y_pert_gauss(:,4))/(2*pi);
semilogy(tspan/T_orb, error_OM, 'b');
grid on
title('\Omega')
ylabel('|\Omega_{cart} - \Omega_{gauss}|/2\pi')
xlabel('Time in numbers of periods')

% pericenter anomaly
subplot(3,2,5)
error_om = abs(PER_AN(:)-Y_pert_gauss(:,5))/(2*pi);
semilogy(tspan/T_orb, error_om, 'b');
grid on
title('\omega')
ylabel('|\omega_{cart} - \omega_{gauss}|/2\pi')
xlabel('Time in numbers of periods')

% true anomaly
subplot(3,2,6)
error_th = abs(rad2deg(TH(:)-Y_pert_gauss(:,6)))./abs(rad2deg(Y_pert_gauss(:,6)));
semilogy(tspan/T_orb, error_th, 'b');
grid on
title('\theta')
ylabel('|\theta_{cart} - \theta_{gauss}|/|\theta_{gauss}|')
xlabel('Time in numbers of periods')


%% Behaviour of keplerian parameters and filtering

figure
subplot(3,2,1)
a_filtered_long = movmean(Y_pert_gauss(:,1), T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,1));
hold on
plot(tspan./T_orb, a_filtered_long, 'k');
grid on                                                                                                   
title('a')

subplot(3,2,2)
e_filtered_long = movmean(Y_pert_gauss(:,2),  T_orb/4);
plot(tspan./T_orb, Y_pert_gauss(:,2));
hold on
plot(tspan./T_orb, e_filtered_long, 'k');
grid on                                                                                                   
title('e')

subplot(3,2,3)
i_filtered_long = movmean(Y_pert_gauss(:,3),  T_orb/4);
plot(tspan./T_orb, Y_pert_gauss(:,3));
hold on
plot(tspan./T_orb, i_filtered_long, 'k');
grid on                                                                                                   
title('i')

subplot(3,2,4)
OM_filtered_long = movmean(Y_pert_gauss(:,4),  T_orb/4);
OM_filtered_sec = movmean(Y_pert_gauss(:,4),  4.5*T_orb);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,4)), 'b');
hold on
plot(tspan./T_orb, rad2deg(OM_filtered_long), 'r');
plot(tspan./T_orb, rad2deg(OM_filtered_sec), 'k--');
grid on                                                                                                   
title('\Omega')

subplot(3,2,5)
om_filtered_long = movmean(Y_pert_gauss(:,5),  T_orb/4);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,5)));
hold on
plot(tspan./T_orb, rad2deg(om_filtered_long), 'k');
grid on                                                                                                   
title('\omega')

subplot(3,2,6)
th_filtered_long = movmean(Y_pert_gauss(:,6),  T_orb/2);
plot(tspan./T_orb, Y_pert_gauss(:,6));
hold on
plot(tspan./T_orb, th_filtered_long, 'k');
grid on    
title('\theta')
hold on
inset2 = axes('Position', [0.1 0.1 0.1 0.1]);
plot(inset2, tspan./T_orb, Y_pert_gauss(:,6));
hold on
plot(inset2, tspan./T_orb, th_filtered_long, 'k')
xlim(inset2, [tspan(100)/T_orb tspan(200)/T_orb]); % Limita l'area di zoom sull'asse x
ylim(inset2, [Y_pert_gauss(100,6) Y_pert_gauss(200,6)]); % Limita l'area di zoom sull'asse y


%%
% Dati di esempio
x = linspace(0, 10, 100);
y = sin(x);

% Crea il plot principale
figure;
plot(x, y);
hold on;

% Definisci l'area di zoom (posizione e dimensioni relative)
zoomPosition = [0.6 0.6 0.25 0.25];
axes('Position', zoomPosition);

% Crea il plot nell'area di zoom
plot(x, y);
xlim([2 4]); % Limita l'area di zoom sull'asse x
ylim([sin(4) sin(2)]); % Limita l'area di zoom sull'asse y

% Aggiungi un riquadro che indica l'area di zoom sul plot principale
rectangle('Position', [2, min(y), 2, max(y)-min(y)], 'EdgeColor', 'r', 'LineWidth', 2);


%% Behaviour as single plots

figure
a_filtered = movmean(Y_pert_gauss(:,1), 10*T_orb);
plot(tspan./T_orb, Y_pert_gauss(:,1), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, a_filtered, 'r');
grid on                                                                                                   
title('a evolution over time')
ylabel('a [km]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered')
%%
figure
e_filtered = movmean(Y_pert_gauss(:,2),  T_orb/5);
e_filtered_secular = movmean(Y_pert_gauss(:,2),  1/0.77e-5);
plot(tspan./T_orb, Y_pert_gauss(:,2), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, e_filtered, 'r');
plot(tspan./T_orb, e_filtered_secular, 'y--', LineWidth=2);
grid on   
title('e evolution over time')
ylabel('e [-]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered', 'Secular behaviour')
inset_e = axes('Position', [0.65 0.15 0.25 0.25]);
plot(inset_e, tspan./T_orb, Y_pert_gauss(:,2), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_e, tspan./T_orb, e_filtered, 'r');
title('Short-period variation filtered')
xlim(inset_e, [20 29]); 
ylim(inset_e, [0.348 0.3495]);
%%
figure
i_filtered = movmean(Y_pert_gauss(:,3),  T_orb/5);
i_filtered_secular = movmean(Y_pert_gauss(:,3),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,3)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(i_filtered), 'r');
plot(tspan./T_orb, rad2deg(i_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('i evolution over time')
ylabel('i [deg]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered', 'Secular behaviour')
inset_i = axes('Position', [0.16 0.61 0.25 0.25]);
plot(inset_i, tspan./T_orb, rad2deg(Y_pert_gauss(:,3)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_i, tspan./T_orb, rad2deg(i_filtered), 'r');
title('Short-period variation filtered')
xlim(inset_i, [60 69]); 
ylim(inset_i, [18.54 18.56]);
%%
figure
OM_filtered = movmean(Y_pert_gauss(:,4),  T_orb/5);
OM_filtered_secular = movmean(Y_pert_gauss(:,4),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,4)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(OM_filtered), 'r');
plot(tspan./T_orb, rad2deg(OM_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('\Omega evolution over time')
ylabel('\Omega [deg]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered', 'Secular behaviour')
inset_OM = axes('Position', [0.2 0.2 0.25 0.25]);
plot(inset_OM, tspan./T_orb, rad2deg(Y_pert_gauss(:,4)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_OM, tspan./T_orb, rad2deg(OM_filtered), 'r');
title('Short-period variation filtered')
xlim(inset_OM, [67 74]); 
ylim(inset_OM, [150 150.1]);
%%
figure
om_filtered = movmean(Y_pert_gauss(:,5),  T_orb/5);
om_filtered_secular = movmean(Y_pert_gauss(:,5),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,5)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(om_filtered), 'r');
plot(tspan./T_orb, rad2deg(om_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('\omega evolution over time')
ylabel('\omega [deg]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered', 'Secular behaviour')
inset_om = axes('Position', [0.15 0.15 0.25 0.25]);
plot(inset_om, tspan./T_orb, rad2deg(Y_pert_gauss(:,5)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_om, tspan./T_orb, rad2deg(om_filtered), 'r');
title('Short-period variation filtered')
xlim(inset_om, [59.5 63]); 
ylim(inset_om, [57.4 57.7]);
%%
figure
th_filtered = movmean(Y_pert_gauss(:,6),  T_orb/5);
th_filtered_secular = movmean(Y_pert_gauss(:,6),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss(:,6)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(th_filtered), 'r');
plot(tspan./T_orb, rad2deg(th_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('\theta evolution over time')
ylabel('\theta [deg]')
xlabel('Time [T]')
xlim([0 1000])
legend('Unfiltered', 'Filtered', 'Secular behaviour')
inset3 = axes('Position', [0.2 0.6 0.25 0.25]);
plot(inset3, tspan./T_orb, rad2deg(Y_pert_gauss(:,6)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset3, tspan./T_orb, rad2deg(th_filtered), 'r');
title('Short-period variation filtered')
xlim(inset3, [59 61]); 
ylim(inset3, 1e+4 *[2.13 2.2]);

% Dati di esempio



%% SIDE NOTES
% - with only J2 effect active we see, as espected from theory:
%      - regression of OMEGA
%      - procession of omega 
% - with J2+3BP effects active we see:
%      - regression of OMEGA
%      - regression of omega 
%      - procession of i









%% 7TH REQUEST
% two objects (non operative satellites) orbiting around Earth
% Data from Space Track (TLE) and NASA Horizon System (ephemeris)

clear all
close all
clc

% Data for Earth as primary planet
mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 12 14 0 0 0]);

%% ENVISAT (Low Earth Orbit)
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


%% IMAGE (High Elliptical Orbit) 
% upload the information from goes13_eph
load('IMAGEeph.mat');

% a_eph_goes = goes12eph(:,10);
% e_eph_goes = goes12eph(:,1);
% i_eph_goes = deg2rad(goes12eph(:,3));
% OM_eph_goes = deg2rad(goes12eph(:,4));
% om_eph_goes = deg2rad(goes12eph(:,5));
% th_eph_goes = deg2rad(goes12eph(:,9));

% a_eph_goes = goes13eph(:,10);
% e_eph_goes = goes13eph(:,1);
% i_eph_goes = deg2rad(goes13eph(:,3));
% OM_eph_goes = deg2rad(goes13eph(:,4));
% om_eph_goes = deg2rad(goes13eph(:,5));
% th_eph_goes = deg2rad(goes13eph(:,9));

a_eph_goes = IMAGEeph(:,10);
e_eph_goes = IMAGEeph(:,1);
i_eph_goes = deg2rad(IMAGEeph(:,3));
OM_eph_goes = deg2rad(IMAGEeph(:,4));
om_eph_goes = deg2rad(IMAGEeph(:,5));
th_eph_goes = deg2rad(IMAGEeph(:,9));


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


