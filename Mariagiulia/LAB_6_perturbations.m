%%LAB_6_PERTURBATIONS

%% J2 perturbation 
clear all
close all
clc

mu_E = astroConstants(13);
R_E = astroConstants(23);
J2 = astroConstants(9);

%day = date2mjd2000([2024, 12, 14, 15, 24, 0)
day = date2mjd2000([2000, 1, 1, 12, 0, 0]);

kep_0 = [7571, 0.01, deg2rad(90), deg2rad(180), deg2rad(180), deg2rad(0)];

[r0, v0] = par2car(kep_0(1), kep_0(2), kep_0(3), kep_0(4), kep_0(5), kep_0(6), mu_E);
T = 2*pi*sqrt(kep_0(1)^3/mu_E);
tspan = linspace(0, 100*T, 100000);
tspan = tspan';

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

[T_ode, S] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(T, k_el, J2, mu_E, R_E, day), mu_E, J2, R_E, day), tspan, kep_0, options);
[T_ODE, Y_cart] = Orbit_Analysis(r0, v0, mu_E, tspan, 'perturbed', 'plot_orbit', 'no_plots', J2, R_E);

[r00,~]  = par2car(S(1, 1), S(1, 2), S(1, 3), S(1, 4), S(1, 5), S(1, 6), mu_E);
p = S(:, 1).*(1-S(:, 2).^2);
r = p./(1+S(:,2).*cos(S(:,6)));

Y_gauss = zeros(size(S, 1), 3);
for i = 1:size(S, 1)
    [rr, v] = par2car(S(i, 1), S(i, 2), S(i, 3), S(i, 4), S(i, 5), S(i, 6), mu_E);
    Y_gauss(i,:) = rr;
end


figure
plot3(Y_gauss(:,1), Y_gauss(:,2), Y_gauss(:,3), 'r');
hold on
plot3(Y_cart(:,1), Y_cart(:,2), Y_cart(:,3), 'b');
axis equal
grid on
legend('gaussian', 'cartesian');


figure
R = Y_cart(:,1:3);
AS = vecnorm(R, 2, 2);
error = r - AS;
plot([1:1:length(error)], error);

error_x = abs(Y_gauss(:,1) - Y_cart(:,1));
error_y = abs(Y_gauss(:,2) - Y_cart(:,2));
error_z = abs(Y_gauss(:,3) - Y_cart(:,3));

figure
plot(tspan, error_x, 'r');
hold on
plot(tspan, error_y, 'b');
plot(tspan, error_z, 'g');

figure
a_filtered = movmean(S(:,1), [10000 10000]);
plot(tspan./T, S(:,1));
hold on
plot(tspan./T, a_filtered, 'k');
grid on                                                                                                   
title('a')

figure
e_filtered = movmean(S(:,2), [10000 10000]);
plot(tspan./T, S(:,2));
hold on
plot(tspan./T, e_filtered, 'k');
grid on
title('e')

figure
i_filtered = movmean(S(:,3), [10000 10000]);
plot(tspan./T, S(:,3));
hold on
plot(tspan./T, i_filtered, 'k');
grid on
axis equal
title('i')

figure
OM_filtered = movmean(S(:,4), [10000 10000]);
plot(tspan./T, rad2deg(S(:,4)));
hold on
%plot(tspan./T, rad2deg(OM_filtered), 'k');
title('\Omega')

figure
om_filtered = movmean(S(:,5), [10000 10000]);
plot(tspan./T, rad2deg(S(:,5)));
hold on
%plot(tspan./T, rad2deg(om_filtered), 'k');
grid on
axis equal
title('\omega')

figure
th_filtered = movmean(S(:,6), [10000 10000]);
plot(tspan./T, S(:,6));
hold on
plot(tspan./T, th_filtered, 'k');
grid on
axis equal
title('\theta')