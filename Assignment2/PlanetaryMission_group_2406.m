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
% DATE: 6/1/2025
%
% -------------------------------------------------------------------------

clc
clear
close all

addpath 'functions'\
addpath 'functions'\timeConversion\time\

%set(0,'defaultfigurecolor',[1 1 1])

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
r_mean_Me   = 5.791e+7;                                 %[km]           Mercury's orbit mean ratio
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
start_date  = date2mjd2000([2041 3 28 0 0 0]);          %[mjd2000]  Starting date of the planetary mission                          

% Orbit initial conditions in cartesian coordinates
[rr_0, vv_0] = par2car (a_0, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_coord_0 = [rr_0; vv_0];

% Repeating GT ratio assigned
m            = 3;                                        %[-]    rotations of the planet
k            = 5;                                        %[-]    revolutions of the satellite

% Options for the integrator
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% number of orbits and points for the simulations
n_orbits = 500;
step = T_orb/500;
tspan = [0:step:n_orbits*T_orb];


%% Propagation and ground track of the unperturbed assigned orbit

% 1 orbit
tspan = [0:step/10:1*T_orb];
% Propagate orbit without perturbation acting on it
[T_1,Y_1] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan, 'non_perturbed');
delta_lambda = rad2deg(T_orb * OMEGA_Me);
fprintf('GT: Delta Lambda = %f \n', delta_lambda);

% Calculate latitude and longitude of the ground track
[~, ~, lon, lat] = groundTrack_2(T_1, OMEGA_Me, 0, 0, Y_1);

% 1 days
tspan_1_day = [0:step/10:24*60*60];
% Propagate orbit without perturbation acting on it
[T_10,Y_10] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan_1_day, 'non_perturbed');
% Calculate latitude and longitude of the ground track
[~, ~, lon_10, lat_10] = groundTrack_2(T_10, OMEGA_Me, 0, 0, Y_10);

% 10 days
tspan_10_days = [0:step:10*24*60*60];
% Propagate orbit without perturbation acting on it
[T_100,Y_100] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan_10_days, 'non_perturbed');
% Calculate latitude and longitude of the ground track
[~, ~, lon_100, lat_100] = groundTrack_2(T_100, OMEGA_Me, 0, 0, Y_100);

% Plot ground track
figure
%subplot(3,1,1)
plotGroundTrack(gca, lon, lat, T_1);
title('Unperturbated Ground Track: 1 period', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,1,2)
figure
plotGroundTrack(gca, lon_10, lat_10, T_10);
title('Unperturbated Ground Track: 1 day', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,1,3)
figure
plotGroundTrack(gca, lon_100, lat_100, T_100);
title('Unperturbated Ground Track: 10 days', Interpreter='latex', FontSize=22, FontWeight='bold');

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
text(1.6*R_Me, 0, 0, 'X', 'FontSize', 9, 'Color', [178/255, 0, 0], Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 0, 1.5*R_Me, 0, 'Color', [178/255, 0, 0], LineWidth=1); 
text(0, 1.6*R_Me,  0, 'Y', 'FontSize', 9, 'Color', [178/255, 0, 0], Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 0, 0, 1.5*R_Me, 'Color',[178/255, 0, 0], LineWidth=1); 
text(0, 0, 1.6*R_Me,  'Z', 'FontSize', 9, 'Color', [178/255, 0, 0], Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 1.5*R_Me*N(1), 1.5*R_Me*N(2), N(3), 'Color', [255/255, 165/255, 0/255], LineWidth=1); 
text(1.5*R_Me*cos(OM_0), 1.5*R_Me*sin(OM_0), 0, 'N', 'FontSize', 9, 'Color', [255/255, 165/255, 0/255], Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 1.5*R_Me*h(1), 1.5*R_Me*h(2), 1.5*R_Me*h(3), 'k', LineWidth=1); 
text( 1.6*R_Me*h(1), 1.6*R_Me*h(2), 1.6*R_Me*h(3), 'h', 'FontSize', 9, 'Color', 'k', Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 1.5*R_Me*e(1), 1.5*R_Me*e(2), 1.5*R_Me*e(3), 'k', LineWidth=1); 
text(1.6*R_Me*e(1), 1.6*R_Me*e(2), 1.6*R_Me*e(3), 'e', 'FontSize', 9, 'Color', 'k', Interpreter='latex', FontSize=12);
quiver3(0, 0, 0, 1.5*R_Me*p(1), 1.5*R_Me*p(2), 1.5*R_Me*p(3), 'k', LineWidth=1); 
text(1.6*R_Me*p(1), 1.6*R_Me*p(2), 1.6*R_Me*p(3), 'p', 'FontSize', 9, 'Color', 'k', Interpreter='latex', FontSize=12);
grid on
axis equal
title('Nominal Orbit', Interpreter='latex', FontSize=14, FontWeight='bold');
legend('Mercury', 'Nominal Orbit' ,'Equatorial Reference Frame', '', '', 'Ascending Node', 'Orbital Reference Frame', '', '', Interpreter='latex', FontSize=12);


%% Compute modified semi-major axis to obtain repeating orbits and groundtrack with a_modif

% modified semi-major axis
n = OMEGA_Me *k/m;
a_modif = (mu_Me/n^2)^(1/3);

% Orbit initial conditions with a_modif in cartesian coordinates
[rr_0_modif, vv_0_modif] = par2car (a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif
T_orb_modif = 2*pi*sqrt(a_modif^3/mu_Me);
delta_lambda_modif = rad2deg(T_orb_modif * OMEGA_Me);

% after k S/C orbits the ground track will repeat equal to itself
tspan_1 = [0:step:k*T_orb_modif];

% Propagate orbit with a_modif without perturbation acting on it
[T_modif,Y_modif] = Orbit_Analysis(rr_0_modif, vv_0_modif, mu_Me, tspan_1, 'non_perturbed');

% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif, lat_modif] = groundTrack_2(T_modif, OMEGA_Me, 0, 0, Y_modif);

% Plot ground track
figure
plotGroundTrack(gca, lon_modif, lat_modif, T_modif);
title('Unperturbed Repeating Ground Track, k=5, m=3', interpreter='latex', FontSize=14);

% Plot Orbit
Mercury_3D
hold on
plotOrbit(Y_1, 2);
plotOrbit(Y_modif, 2);
r_SOI = r_mean_Me * (mu_Me/mu_Sun)^(2/5); 
[rr_A_modif, ~] = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);
plot3(rr_A_modif(1), rr_A_modif(2), rr_A_modif(3), 'k*', LineWidth=1)
plotSOI(r_SOI)
grid on
axis equal
title('Unperturbed orbit for the repeating ground track', Interpreter='latex', FontSize=14, FontWeight='bold');
legend('Mercury', 'Nominal Orbit', 'Modified orbit for RGT', 'Apocenter','SOI', Interpreter='latex', FontSize=12)


%% Check: is the apocenter of the modified orbit inside the sphere of interest of Mercury?

% apocenter and pericenter position vector
[rr_A_modif, ~] = par2car (a_modif, e_0, i_0, OM_0, om_0, pi, mu_Me);
r_apo = norm(rr_A_modif);
[rr_P_modif, ~] = par2car (a_modif, e_0, i_0, OM_0, om_0, 0, mu_Me);
r_peri = norm(rr_P_modif);

r_SOI = r_mean_Me * (mu_Me/mu_Sun)^(2/5);        %sphere of influence of Mercury
n_perc = 0.5;                                    % limit on the percentage of SOI for apocenter

a_rep = zeros(30,30);
RATIO = zeros(30,30);

flag = 1;
if r_apo > r_SOI 
    fprintf('Apocenter radius is larger than the sphere of interest of Mercury\n');
    for K = [1:1:30]                                                        % span for k
        for M = [1:1:30]                                                    % span for m
            k_m_ratio = K/M;                                                % k/m ration evaluated
            a_modif_new = (mu_Me*(M/(K*OMEGA_Me))^2)^(1/3);                 % calculation of the semi-major axis with k/m
            a_rep(K,M) = a_modif_new;
            RATIO(K,M) = k_m_ratio;
            ra_new = a_modif_new * (1+e_0);
            %if a_modif_new <= ((n_perc*r_SOI)/(1+e_0))                      % verify if a is below a certain threshold that  
            if ra_new <= (n_perc*r_SOI)       
                a_rep_valid(flag) = a_modif_new;                            % depends on the constraint of r_apo to be below the 
                RATIO_valid(flag) = k_m_ratio;                              % 50% of the SOI                
                k_valid(flag) = K;
                m_valid(flag) = M;
                flag = flag+1;   
            end
        end
    end
else 
    fprintf('Apocenter radius is smaller than the sphere of interest of Mercury\n');
end

% first one has the lowest k and it has been chosen
a_modif_new = a_rep_valid(1);
[r_apo_new, ~] = par2car(a_modif_new, e_0, i_0, OM_0, om_0, pi, mu_Me);
ra_new_new = a_modif_new * (1+e_0);
r_apo_new_norm = norm(r_apo_new);
perc = r_apo_new_norm/r_SOI;
fprintf('For the new repeating GT, the condition of r_apo <= %.2f *r_SOI is respected: r_apo_new/r_SOI = %f \n', n_perc, perc);
k_new = k_valid(1);
m_modif_new = m_valid(1);
fprintf('Valid value of semi-major axix: %.3f \n', a_modif_new);
fprintf('Valid value of k/m: %.1f, where k = %.1f and m = %.1f \n', k_new/m_modif_new, k_new, m_modif_new);

% Orbit initial conditions with a_modif_new in cartesian coordinates
[rr_0_modif_new, vv_0_modif_new] = par2car (a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);

% Orbit period with a_modif_new
T_orb_modif_new  = 2*pi*sqrt(a_modif_new^3/mu_Me);
fprintf('Valid value of period: %.2f \n', T_orb_modif_new);
delta_lambda_new = rad2deg(T_orb_modif_new * OMEGA_Me);
fprintf('New Repating GT: Delta Lambda = %.2f \n\n', delta_lambda_new);

% 1 orbit
tspan_2 = [0:step:T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
[T_modif_new,Y_modif_new] = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_2, 'non_perturbed');
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_new, lat_modif_new] = groundTrack_2(T_modif_new, OMEGA_Me, 0, 0, Y_modif_new);

% 15 orbits (k_new = 15)
tspan_15 = [0:step:k_new*T_orb_modif_new];
% Propagate orbit with a_modif without perturbation acting on it
[T_modif_new_15,Y_modif_new_15] = Orbit_Analysis(rr_0_modif_new, vv_0_modif_new, mu_Me, tspan_15, 'non_perturbed');
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_new_15, lat_modif_new_15] = groundTrack_2(T_modif_new_15, OMEGA_Me, 0, 0, Y_modif_new_15);

% Plot ground track
figure
%subplot(2,1,1)
plotGroundTrack(gca, lon_modif_new, lat_modif_new, T_modif_new);
title('Unperturbed Repeating Ground Track, k=15, m=1: 1 orbit', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(2,1,2)
figure
plotGroundTrack(gca, lon_modif_new_15, lat_modif_new_15, T_modif_new_15);
title('Unperturbed Repeating Ground Track, k=15, m=1: 15 orbits', Interpreter='latex', FontSize=22, FontWeight='bold');

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
title('Unperturbed orbit for the repeating ground track and the new one', Interpreter='latex', FontSize=14, FontWeight='bold');
legend('Mercury', 'Nominal Orbit','Repeating orbit', 'Apocenter', 'New repeating orbit', 'SOI', Interpreter='latex', FontSize=12);

% Orbit plot: nominal and new repeating
Mercury_3D
hold on
plotOrbit(Y_1, 2)
plotOrbit(Y_modif_new, 2)
grid on
axis equal
title('Unperturbed orbit for the new repeating ground track', Interpreter='latex', FontSize=14, FontWeight='bold');
legend('Mercury', 'Nominal Orbit', 'New repeating orbit', Interpreter='latex', FontSize=12);


%% Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator, Cartesian
% The following section is here to demonstrate that outside Mercury's SOI
% the third-body perturbation model can be applied
                                                                               % SIDE NOTES: T_orb_modif is 105 times larger than nominal,
T_orb_modif  = 2*pi*sqrt(a_modif^3/mu_Me);                                     % solar influence gets too big 
time_check = T_orb_modif/T_orb;                                                
tspan_1 = [0:step:50*T_orb_modif];                                             % it encounter some sort of discontinuity 

% cartesia initial condition
[rr_0_test, vv_0_test] = par2car(a_modif, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
c_0_modif = [rr_0_test; vv_0_test];

% cartesian propagator
[T_pert_cart_modif, Y_pert_cart_modif] = ode113(@(t, k_el) eq_motion_cartesian(t, k_el, @(t, k_el) acc_pert_function_cartesian(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, c_0_modif, options);

% Orbit Plot 
Mercury_3D
hold on
plot3(Y_pert_cart_modif(:,1), Y_pert_cart_modif(:,2), Y_pert_cart_modif(:,3))
scatterOrbit(Y_pert_cart_modif, tspan_1)
grid on
axis equal
title('orbit - cartesian propagator', Interpreter='latex', FontSize=14, FontWeight='bold');


%% CORRECTED Reapeating GT with perturbations: J2 and 3BP - Sun influence: Gaussian propagator

% new time span
tspan_1 = [0:step:k_new*T_orb_modif_new];

% set of initial conditions with the new semi-major axis
k_el_0_modif = [a_modif_new; e_0; i_0; OM_0; om_0; theta_0];
[r_mod_0, v_mod_0] = par2car(a_modif_new, e_0, i_0, OM_0, om_0, theta_0, mu_Me);
cart_0_modif = [r_mod_0; v_mod_0];

% propagation with Gauss Planetary equations
[T_pert_gauss_modif, Y_pert_gauss_modif] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, k_el_0_modif, options);

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
title('Ground Track of the perturbated orbit ($a_{modif}$) - J2 and 3BP', Interpreter='latex', FontSize=14, FontWeight='bold');
ylabel(colorbar,'Time');

%Orbit Plot
Mercury_3D
hold on
scatterOrbit(R_pert_gauss_modif, tspan_1);
grid on
axis equal
title('orbit - gaussian propagator', interpreter='latex', FontSize=14, FontWeight='bold');

% Propagation with Cartesian Equation of Motion 
% Computational time: higher than cartesian, but it does not require the
% conversion from keplerian to cartesian for the ground track plot

% 1 orbit
tspan_1 = [0:step:T_orb_modif_new];
[T_pert_cart_modif, Y_pert_cart_modif] = ode113(@(t, cc) eq_motion_cartesian(t, cc, @(t, cc) acc_pert_function_cartesian(t, cc, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_1, cart_0_modif, options);
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pertc, lat_modif_pertc] = groundTrack_2(T_pert_cart_modif, OMEGA_Me, 0, 0, Y_pert_cart_modif);

% 15 orbits
tspan_15 = [0:step:k_new*T_orb_modif_new];
[T_pert_cart_modif_15, Y_pert_cart_modif_15] = ode113(@(t, cc) eq_motion_cartesian(t, cc, @(t, cc) acc_pert_function_cartesian(t, cc, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_15, cart_0_modif, options);
% Calculate latitude and longitude of the ground track with a_modif
[~, ~, lon_modif_pertc_15, lat_modif_pertc_15]    = groundTrack_2(T_pert_cart_modif_15, OMEGA_Me, 0, 0, Y_pert_cart_modif_15);

% Plot ground track perturbed and unperturbed for comparison, as single figures or subplot
figure
%subplot(2,2,1)
plotGroundTrack(gca, lon_modif_new, lat_modif_new, T_modif_new);
title('Unperturbed Repeating Ground Track, k=15, m=1: 1 orbit', interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(2,2,3)
figure
plotGroundTrack(gca, lon_modif_new_15, lat_modif_new_15, T_modif_new_15);
title('Unperturbed Repeating Ground Track, k=15, m=1: 15 orbits', interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(2,2,2)
figure
plotGroundTrack(gca, lon_modif_pertc, lat_modif_pertc, T_pert_cart_modif);
title('Perturbed Repeating Ground Track, k=15, m=1: 1 orbit', interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(2,2,4)
figure
plotGroundTrack(gca, lon_modif_pertc_15, lat_modif_pertc_15, T_pert_cart_modif_15);
title('Perturbed Repeating Ground Track, k=15, m=1: 15 orbits', interpreter='latex', FontSize=22, FontWeight='bold');


%% Introduce perturbations: J2 and 3BP - Sun influence: Gaussian propagator and Cartesian Propagator

% Nominal Initial Conditions
k_el_0      = [a_0;e_0;i_0;OM_0;om_0;theta_0];  

% new tspan appropriate to see the perturbations for the 3D plot
tspan = [0:step:1000*T_orb];

% Gaussian Propagator
tic
[T_pert_gauss, Y_pert_gauss] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, k_el_0, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation for 1000 orbits: %.4f s. \n', elapsedtime_gauss)

% matrix of the position vector 
R_pert_gauss = zeros(size(Y_pert_gauss,1), 3);
for j = 1:size(Y_pert_gauss, 1)
    [rr, ~] = par2car(Y_pert_gauss(j, 1), Y_pert_gauss(j, 2), Y_pert_gauss(j, 3), Y_pert_gauss(j, 4), Y_pert_gauss(j, 5), Y_pert_gauss(j, 6), mu_Me);
    rr = rr';
    R_pert_gauss(j,:) = rr;
end

% Plot Orbit - Gauss
Mercury_3D
hold on
scatterOrbit(R_pert_gauss, tspan);
grid on
axis equal
title('Perturbed orbit over 1000 periods', Interpreter='latex', FontSize=14, FontWeight='bold');
legend('Mercury', 'Perturbed orbit', Interpreter='latex', FontSize=12)

% Cartesian Propagator
tic
[T_pert_cart, Y_pert_cart] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation over 1000 orbits: %.4f s. \n', elapsedtime_cart);

% Plot Orbit - Cartesian
Mercury_3D
hold on
scatterOrbit(Y_pert_cart, tspan);
grid on
axis equal
title('orbit revolution over 1000 periods - cartesian propagator', Interpreter='latex', FontSize=14, FontWeight='bold');


%% Plot of perturbed orbit ground track from Gaussian and Cartesian propagators

% appropriate number of orbits to analize keplerian parameters behaviour
% for the ground track
n_orb_GT = 100;

% more appropriate tspan for the GT
tspan_GT = [0:step:n_orb_GT*T_orb];

[T_pert_gauss_GT, Y_pert_gauss_GT] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT, k_el_0, options);

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
title('Ground Track of the perturbated orbit (Gauss) for 5 orbits: J2 and 3BP-Sun', Interpreter='latex', FontSize=14, FontWeight='bold');


%% Cartesian Propagator
% even if sligly slower than cartesian for this time span, it eneables to
% allocate less memory to save both keplerian parameters and posizion and
% velocity due to the conversion necessary for the plot of the ground track

% 1 orbit
tspan_GT = [0:step/10:T_orb];
tic
[T_pert_cart_GT, Y_pert_cart_GT] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation GT 1 period: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT, lat_pert_cart_GT] = groundTrack_2(T_pert_cart_GT, OMEGA_Me, 0, 0, Y_pert_cart_GT);

% 1 day
tspan_GT_1d = [0:step/10:24*60*60];
tic
[T_pert_cart_GT_1d, Y_pert_cart_GT_1d] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT_1d, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation GT 1 day: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT_1d, lat_pert_cart_GT_1d] = groundTrack_2(T_pert_cart_GT_1d, OMEGA_Me, 0, 0, Y_pert_cart_GT_1d);

% % 10 days
% tspan_10_days = [0:step:10*24*60*60];
% % Propagate orbit without perturbation acting on it
% [T_100,Y_100] = Orbit_Analysis(rr_0, vv_0, mu_Me, tspan_10_days, 'non_perturbed');
% % Calculate latitude and longitude of the ground track
% [~, ~, lon_100, lat_100] = groundTrack_2(T_100, OMEGA_Me, 0, 0, Y_100);

% 10 days
tspan_GT_10d = [0:step:10*24*60*60];
tic
[T_pert_cart_GT_10d, Y_pert_cart_GT_10d] = ode113(@(t, cart_coord) eq_motion_cartesian(t, cart_coord, @(t, cart_coord) acc_pert_function_cartesian(t, cart_coord, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan_GT_10d, cart_coord_0, options);
elapsedtime_cart = toc;
fprintf ('Cartesian simulation GT 10 days: %.4f s. \n', elapsedtime_cart);
% gt
[~, ~, lon_pert_cart_GT_10d, lat_pert_cart_GT_10d] = groundTrack_2(T_pert_cart_GT_10d, OMEGA_Me, 0, 0, Y_pert_cart_GT_10d);

% Plot ground track perturbed and unperturbed for comparison, as single figures or subplots
figure
%subplot(3,2,1)
plotGroundTrack(gca, lon, lat, T_1);
title('Unperturbated Ground Track: 1 period', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,2,3)
figure
plotGroundTrack(gca, lon_10, lat_10, T_10);
title('Unperturbated Ground Track: 1 day', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,2,5)
figure
plotGroundTrack(gca, lon_100, lat_100, T_100);
title('Unperturbated Ground Track: 10 days', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,2,2)
figure
plotGroundTrack(gca,lon_pert_cart_GT, lat_pert_cart_GT, T_pert_cart_GT);
title('Perturbed Ground Track: 1 period', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,2,4)
figure
plotGroundTrack(gca,lon_pert_cart_GT_1d, lat_pert_cart_GT_1d, T_pert_cart_GT_1d);
title('Perturbed Ground Track: 1 day', Interpreter='latex', FontSize=22, FontWeight='bold');
%subplot(3,2,6)
figure
plotGroundTrack(gca,lon_pert_cart_GT_10d, lat_pert_cart_GT_10d, T_pert_cart_GT_10d);
title('Perturbed Ground Track: 10 days', Interpreter='latex', FontSize=22, FontWeight='bold');


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

% unwrap of the vectors of angles
RAAN = unwrap(RAAN);
PER_AN = unwrap(PER_AN);
TH = unwrap(TH);

% semi-major axis 
figure
subplot(3,2,1)
error_a = abs(A(:)-Y_pert_gauss(:,1))/a_0;
a_e_max = max(error_a);
semilogy(tspan/T_orb, error_a, 'b');
grid on
xlim([0 1000])
title('a', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|a_{cart}$ - $a_{gauss}|$/$a_0$', Interpreter='latex', FontSize=12)
xlabel('Time [T]', Interpreter='latex', FontSize=12)

% eccentricity 
subplot(3,2,2)
error_e = abs(E(:)-Y_pert_gauss(:,2));
e_e_max = max(error_e);
semilogy(tspan/T_orb, error_e, 'b');
grid on
xlim([0 1000])
title('e', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|e_{cart}$ - $e_{gauss}|$', Interpreter='latex', FontSize=12)
xlabel('Time [T]', Interpreter='latex', FontSize=12)

% inclination  
subplot(3,2,3)
error_i = abs(I(:)-Y_pert_gauss(:,3));
i_e_max = max(error_i);
semilogy(tspan/T_orb, error_i, 'b');
grid on
xlim([0 1000])
title('i', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|i_{cart}$ - $i_{gauss}|$', Interpreter='latex', FontSize=12)
xlabel('Time [T]', Interpreter='latex', FontSize=12)

% RAAN  
subplot(3,2,4)
error_OM = abs(RAAN(:)-Y_pert_gauss(:,4))/(2*pi);
OM_e_max = max(error_OM);
semilogy(tspan/T_orb, error_OM, 'b');
grid on
xlim([0 1000])
title('$\Omega$', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|\Omega_{cart}$ - $\Omega_{gauss}|$/2$\pi$', Interpreter='latex', FontSize=12)
xlabel('Time [T]', Interpreter='latex', FontSize=12)

% pericenter anomaly
subplot(3,2,5)
error_om = abs(PER_AN(:)-Y_pert_gauss(:,5))/(2*pi);
om_e_max = max(error_om);
semilogy(tspan/T_orb, error_om, 'b');
grid on
xlim([0 1000])
title('$\omega$', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|\omega_{cart}$ - $\omega_{gauss}|$/2$\pi$', Interpreter = 'latex')
xlabel('Time [T]', Interpreter='latex', FontSize=12)

% true anomaly
subplot(3,2,6)
error_th = abs(rad2deg(TH(:)-Y_pert_gauss(:,6)))./abs(rad2deg(Y_pert_gauss(:,6)));
th_e_max = max(error_th);
semilogy(tspan/T_orb, error_th, 'b');
grid on
xlim([0 1000])
title('$\theta$', Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$|\theta_{cart}$ - $\theta_{gauss}|$/$|\theta_{gauss}|$', Interpreter='latex', FontSize=12)
xlabel('Time [T]', Interpreter='latex', FontSize=12)

error = [a_e_max, e_e_max, i_e_max, OM_e_max, om_e_max, th_e_max];

%% Behaviour of keplerian parameters and filtering

% new tspan appropriate to see the effect of perturbation
tspan = [0:step:1000*T_orb];
%tspan = [0:step:20000*T_orb];

% Gaussian Propagator
[T_pert_gauss, Y_pert_gauss_1000] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_Me, mu_Me, R_Me, start_date, eps_Me, 1), mu_Me), tspan, k_el_0, options);

figure
a_filtered = movmean(Y_pert_gauss_1000(:,1), 20*T_orb);
plot(tspan./T_orb, Y_pert_gauss_1000(:,1), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, a_filtered, 'r');
grid on                                                                                                   
title('Semi-major axis evolution over time',  Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('a [km]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1',  interpreter='latex', fontsize=12)

figure
e_filtered = movmean(Y_pert_gauss_1000(:,2),  T_orb/5);
e_filtered_secular = movmean(Y_pert_gauss_1000(:,2),  1/0.77e-5);
plot(tspan./T_orb, Y_pert_gauss_1000(:,2), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, e_filtered, 'r');
plot(tspan./T_orb, e_filtered_secular, 'y--', LineWidth=2);
grid on   
title('Eccentricity evolution over time',  Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('e [-]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1', 'Filter 2',  interpreter='latex', fontsize=12)
inset_e = axes('Position', [0.65 0.15 0.25 0.25]);
plot(inset_e, tspan./T_orb, Y_pert_gauss_1000(:,2), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_e, tspan./T_orb, e_filtered, 'r');
%title('Short-period variation',  Interpreter='latex', FontSize=14, FontWeight='bold')
grid on
xlim(inset_e, [20 29]); 
ylim(inset_e, [0.348 0.3495])

figure
i_filtered = movmean(Y_pert_gauss_1000(:,3),  T_orb/5);
i_filtered_secular = movmean(Y_pert_gauss_1000(:,3),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,3)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(i_filtered), 'r');
plot(tspan./T_orb, rad2deg(i_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('Inclination evolution over time',  interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('i [deg]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1', 'Filter 2',  interpreter='latex', FontSize=14)
inset_i = axes('Position', [0.16 0.63 0.25 0.25]);
plot(inset_i, tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,3)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_i, tspan./T_orb, rad2deg(i_filtered), 'r');
%title('Short-period variation',  interpreter='latex')
grid on
xlim(inset_i, [60 69]); 
ylim(inset_i, [18.54 18.56]);

figure
OM_filtered = movmean(Y_pert_gauss_1000(:,4),  T_orb/5);
OM_filtered_secular = movmean(Y_pert_gauss_1000(:,4),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,4)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(OM_filtered), 'r');
plot(tspan./T_orb, rad2deg(OM_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('RAAN evolution over time',  Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$\Omega$ [deg]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1', 'Filter 2', interpreter='latex', fontsize=12)
inset_OM = axes('Position', [0.2 0.2 0.25 0.25]);
plot(inset_OM, tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,4)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_OM, tspan./T_orb, rad2deg(OM_filtered), 'r');
%title('Short-period variation',  interpreter='latex')
grid on
xlim(inset_OM, [67 74]); 
ylim(inset_OM, [150 150.1]);

figure
om_filtered = movmean(Y_pert_gauss_1000(:,5),  T_orb/5);
om_filtered_secular = movmean(Y_pert_gauss_1000(:,5),  1/0.77e-5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,5)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(om_filtered), 'r');
plot(tspan./T_orb, rad2deg(om_filtered_secular), 'y--', LineWidth=2);
grid on                                                                                                   
title('Argument of pericenter evolution over time',  Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$\omega$ [deg]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1', 'Filter 2',  interpreter='latex', FontSize=12)
inset_om = axes('Position', [0.63 0.6 0.25 0.25]);
plot(inset_om, tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,5)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset_om, tspan./T_orb, rad2deg(om_filtered), 'r');
%title('Short-period variation',  interpreter='latex')
grid on
xlim(inset_om, [59.5 63]); 
ylim(inset_om, [57.4 57.7]);

figure
th_filtered = movmean(Y_pert_gauss_1000(:,6),  T_orb/5);
plot(tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,6)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(tspan./T_orb, rad2deg(th_filtered), 'r');
grid on                                                                                                   
title('True anomaly evolution over time',  Interpreter='latex', FontSize=14, FontWeight='bold')
ylabel('$\theta$ [deg]',  interpreter='latex', FontSize=14)
xlabel('Time [T]',  interpreter='latex', FontSize=14)
xlim([0 1000])
legend('Unfiltered', 'Filter 1',  interpreter='latex', FontSize=12)
inset3 = axes('Position', [0.2 0.6 0.25 0.25]);
plot(inset3, tspan./T_orb, rad2deg(Y_pert_gauss_1000(:,6)), 'Color', [1, 0.6, 0], LineWidth=1);
hold on
plot(inset3, tspan./T_orb, rad2deg(th_filtered), 'r');
%title('Short-period variation',  interpreter='latex')
grid on
xlim(inset3, [59 61]); 
ylim(inset3, 1e+4 *[2.13 2.2]);


%% REAL DATA ANALYSIS
% two objects (non operative satellites) orbiting around Earth
% Data from Space Track (TLE) and NASA Horizon System (ephemerides)

%% IMAGE (High Elliptical Orbit) 
% upload the information from IMAGE1yeareph by NASA/JPL Horizon
% it's a matrix of 12 columns containing all the info from the ephemeris
% generator

load('IMAGE1yeareph.mat');

a_eph_image = IMAGE1yeareph(:,10);
e_eph_image = IMAGE1yeareph(:,1);
i_eph_image = deg2rad(IMAGE1yeareph(:,3));
OM_eph_image = deg2rad(IMAGE1yeareph(:,4));
om_eph_image = deg2rad(IMAGE1yeareph(:,5));
th_eph_image = deg2rad(IMAGE1yeareph(:,9));

% Behaviour
tspan3 = [0:60*24:365*24*60*60];

% unwrap the angles
OM_eph_image = unwrap(OM_eph_image);
om_eph_image = unwrap(om_eph_image);
th_eph_image = unwrap(th_eph_image);

% to have the behaviour expressed in days
vect = tspan3./(60*60*24);

% plot of the behavior of the keplerian elements to compare
figure
subplot(3,2,1)
plot(vect, a_eph_image(:));
grid on
xlim([0 365])
title('Semi-major axis evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('a (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('a [km]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex',  FontSize=14) 

subplot(3,2,2)
plot(vect, e_eph_image(:));
grid on
xlim([0 365])
title('Eccentricity evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('e (ephemerides)', interpreter='latex',FontSize=14) 
ylabel('e [-]', interpreter='latex',FontSize=14) 
xlabel('Time [day]', interpreter='latex',FontSize=14) 

subplot(3,2,3)
plot(vect, rad2deg(i_eph_image(:)));
grid on
xlim([0 365])
title('Inclination evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('i (ephemerides)', interpreter='latex',FontSize=14) 
ylabel('i [deg]', interpreter='latex',FontSize=14) 
xlabel('Time [day]', interpreter='latex',FontSize=14) 

subplot(3,2,4)
plot(vect, rad2deg(OM_eph_image(:)));
grid on
xlim([0 365])
title('RAAN evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\Omega$ (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('$\Omega$ [deg]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 

subplot(3,2,5)
plot(vect, rad2deg(om_eph_image(:)));
grid on
xlim([0 365])
title('Argument of pericenter evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\omega$ (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('$\omega$ [deg]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 


subplot(3,2,6)
plot(vect,  rad2deg(th_eph_image(:)));
grid on
xlim([0 365])
title('True anomaly evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\theta$ (ephemerides)', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 
ylabel('$\theta$ [deg]', interpreter='latex', FontSize=14) 


%%  Alouette-2 (Low Earth Orbit)
% upload the information from Alo21year by NASA/JPL Horizon
% load the matrix 
% it's a matrix of 12 columns containing all the info from the ephemeris generator

load('Alo21year1.mat');

a_eph_Alo2 = Alo21year1(:,10);
e_eph_Alo2 = Alo21year1(:,1);
i_eph_Alo2 = deg2rad(Alo21year1(:,3));
OM_eph_Alo2 = deg2rad(Alo21year1(:,4));
om_eph_Alo2 = deg2rad(Alo21year1(:,5));
th_eph_Alo2 = deg2rad(Alo21year1(:,9));

tspan3 = [0:60*24:365*24*60*60];

% Behaviour

% unwrap the angles
OM_eph_Alo2 = unwrap(OM_eph_Alo2);
om_eph_Alo2 = unwrap(om_eph_Alo2);
th_eph_Alo2 = unwrap(th_eph_Alo2);

% to have the behaviour expressed in days
vect = tspan3./(60*60*24);

% plot of the behavior of the keplerian elements to compare
figure
subplot(3,2,1)
plot(vect, a_eph_Alo2(:));
grid on
xlim([0 365])
ylim([7915 7960])
title('Semi-major axis evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('a (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('a [km]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex',  FontSize=14) 

subplot(3,2,2)
plot(vect, e_eph_Alo2(:));
grid on
xlim([0 365])
ylim([0.125 0.14])
title('Eccentricity evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('e (ephemerides)', interpreter='latex',FontSize=14) 
ylabel('e [-]', interpreter='latex',FontSize=14) 
xlabel('Time [day]', interpreter='latex',FontSize=14) 

subplot(3,2,3)
plot(vect, rad2deg(i_eph_Alo2(:)));
grid on
xlim([0 365])
title('Inclination evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('i (ephemerides)', interpreter='latex',FontSize=14) 
ylabel('i [deg]', interpreter='latex',FontSize=14) 
xlabel('Time [day]', interpreter='latex',FontSize=14) 

subplot(3,2,4)
plot(vect, rad2deg(OM_eph_Alo2(:)));
grid on
xlim([0 365])
title('RAAN evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\Omega$ (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('$\Omega$ [deg]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 

subplot(3,2,5)
plot(vect, rad2deg(om_eph_Alo2(:)));
grid on
xlim([0 365])
title('Argument of pericenter evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\omega$ (ephemerides)', interpreter='latex', FontSize=14) 
ylabel('$\omega$ [deg]', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 

subplot(3,2,6)
plot(vect,  rad2deg(th_eph_Alo2(:)));
grid on
xlim([0 365])
title('True anomaly evolution over time', interpreter='latex', FontSize=14, FontWeight='bold')
legend('$\theta$ (ephemerides)', interpreter='latex', FontSize=14) 
xlabel('Time [day]', interpreter='latex', FontSize=14) 
ylabel('$\theta$ [deg]', interpreter='latex', FontSize=14) 


