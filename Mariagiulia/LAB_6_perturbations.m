%%LAB_6_PERTURBATIONS

%% J2 perturbation 
clear all
close all
clc

mu_E = astroConstants(13);
R_E = astroConstants(23);
J2 = astroConstants(9);

kep_0 = [7571, 0.01, deg2rad(87.9), deg2rad(180), deg2rad(180), deg2rad(0)];
T = 2*pi*sqrt(kep_0(1)^3/mu_E);
tspan = linspace(0, 100*T, 100000);
tspan = tspan';

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

[T_ode, S] = ode113(@(t, k_el) eq_motion_Gauss_SRW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2, mu_E, R_E), mu_E, J2, R_E), tspan, kep_0, options);

p = S(:, 1).*(1-S(:, 2));
r = p./(1-S(:,2).*cos(S(:,6)));

radius = zeros(size(S,1), 3);

for i = 1:size(S, 1)
    [rr, v] = par2car(S(i, 1), S(i, 2), S(i, 3), S(i, 4), S(i, 5), S(i, 6), mu_E);
    radius(i,:) = rr;
end

plot3(radius(:,1), radius(:,2), radius(:,3));

% %%
% plot3(r(:,1), r(:,2), r(:,3))