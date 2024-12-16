%%LAB_6_PERTURBATIONS

clear all
close all
clc

mu_E = astroConstants(13);
R_E = astroConstants(23);
J2 = astroConstants(9);
date = 0;

kep_0 = [7571, 0.01, deg2rad(87.9), deg2rad(180), deg2rad(180), deg2rad(0)];
T_orb = 2*pi*sqrt(kep_0(1)^3/mu_E);
tspan = linspace(0, 10000*T_orb, 100000);
tspan = tspan';

%% Keplerian propagation
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
[T_ode, S] = ode113(@(t, k_el) eq_motion_Gauss_SRW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2, mu_E, R_E,date), mu_E, J2, R_E,date), tspan, kep_0, options);

p = S(:, 1).*(1-S(:, 2).^2);
r = p./(1-S(:,2).*cos(S(:,6)));

Y_kep = zeros(size(S,1), 3);

for i = 1:size(S, 1)
    [rr, v] = par2car(S(i, 1), S(i, 2), S(i, 3), S(i, 4), S(i, 5), S(i, 6), mu_E);
    Y_kep(i,:) = rr;
end

%% Cartesian propagation

% Cartesian propagation
[r0,v0] = par2car(7571, 0.01, deg2rad(87.9), deg2rad(180), deg2rad(180), deg2rad(0),mu_E);
[ T, Y_car ] = orbit_propagator_J2(r0 , v0 , mu_E , 2 , options, 100000 , 10000, J2 , R_E );

%% Plots

% Plot propagated orbits
figure
plot3(Y_kep(:,1), Y_kep(:,2), Y_kep(:,3));
hold on
plot3(Y_car(:,1), Y_car(:,2), Y_car(:,3));
grid on
axis equal
title('Orbit representation');

% Plot errors between the two propagation
error1 = abs(Y_kep(:,1)-Y_car(:,1));
error2 = abs(Y_kep(:,2)-Y_car(:,2));
error3 = abs(Y_kep(:,3)-Y_car(:,3));

figure
semilogy(tspan,error1);
hold on
semilogy(tspan,error2);
hold on
semilogy(tspan,error3);
grid on
title('orbit propagation errors');




