%% TEST PERTURBATION J2
% this scripts tests the propagators
% - orbit_propagator_J2 gives bigger errors

%% J2 perturbation 
clear all
close all
clc

mu_E = astroConstants(13);
R_E = astroConstants(23);
J2 = astroConstants(9);

%day = date2mjd2000([2024, 12, 14, 15, 24, 0)
day = date2mjd2000([2000, 1, 1, 12, 0, 0]);

kep_0 = [7571, 0.01, deg2rad(87.9), deg2rad(180), deg2rad(180), deg2rad(0)];

[r0, v0] = par2car(kep_0(1), kep_0(2), kep_0(3), kep_0(4), kep_0(5), kep_0(6), mu_E);
T = 2*pi*sqrt(kep_0(1)^3/mu_E);
tspan = [0:1:100*T];

%%
coord_0 = [r0; v0];
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 
tilt = astroConstants(63);

[T_ode, S] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2, mu_E, R_E, day, tilt), mu_E), tspan, kep_0, options);
%[T_ode2, Y] = ode113(@(t2, coord) eq_motion_cartesian(t2, coord, @(t2, coord) acc_pert_function_cartesian(t2, coord, J2, mu_E, R_E, day, tilt), mu_E), tspan, coord_0, options);
%[T_ode2,Y]                = orbit_propagator_J2(r0, v0 , mu_E, 2 , options, 655603, 100, J2, R_E);
[T_ode2, Y] = Orbit_Analysis(r0, v0, mu_E, tspan, 'perturbed', J2, R_E);

for i = 1:size(S,1)
    while S(i,6)>2*pi
        S(i,6) = S(i,6) - 2*pi;
    end
end

for i = 1:size(S, 1)
    [a_c, e_c, i_c, OM_c, om_c, th_c] = car2par(Y(i,1:3), Y(i,4:6), mu_E);
    A(i) = a_c;
    E(i) = e_c;
    I(i) = i_c;
    RAAN(i) = OM_c;
    PER_AN(i) = om_c;
    TH(i) = th_c;
end

%% ERRORS

% semi-major axis 
figure
error_a = abs(A(:)-S(:,1))/7571;
semilogy(tspan/T, error_a, 'b');
grid on
title('a')

% eccentricity 
figure
error_e = abs(E(:)-S(:,2));
semilogy(tspan/T, error_e, 'b');
grid on
title('e')

% inclination  
figure
error_i = abs(I(:)-S(:,3));
semilogy(tspan/T, error_i, 'b');
grid on
title('i')

% RAAN  
figure
error_OM = abs(RAAN(:)-S(:,4))/(2*pi);
semilogy(tspan/T, error_OM, 'b');
grid on
title('\Omega')

% pericenter anomaly
figure
error_om = abs(PER_AN(:)-S(:,5))/(2*pi);
semilogy(tspan/T, error_om, 'b');
grid on
title('\omega')

% true anomaly
figure
error_th = abs(rad2deg(TH(:)-S(:,6)))./rad2deg(S(:,6));
semilogy(tspan/T, error_th, 'b');
grid on
title('\theta')

%%
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

%%

figure
R = Y_cart(:,1:3);
AS = vecnorm(R, 2, 2);
error = r - AS;
plot([1:1:length(error)], error);
%%
error_x = abs(Y_gauss(:,1) - Y_cart(:,1));
error_y = abs(Y_gauss(:,2) - Y_cart(:,2));
error_z = abs(Y_gauss(:,3) - Y_cart(:,3));

figure
plot(tspan, error_x, 'r');
hold on
plot(tspan, error_y, 'b');
plot(tspan, error_z, 'g');

%%
figure
a_filtered = movmean(S(:,1), [10000 10000]);
plot(tspan./T, S(:,1));
hold on
plot(tspan./T, a_filtered, 'k');

figure
e_filtered = movmean(S(:,2), [10000 10000]);
plot(tspan./T, S(:,2));
hold on
plot(tspan./T, e_filtered, 'k');

figure
i_filtered = movmean(S(:,3), [10000 10000]);
plot(tspan./T, S(:,3));
hold on
plot(tspan./T, i_filtered, 'k');

figure
OM_filtered = movmean(S(:,4), [10000 10000]);
plot(tspan./T, rad2deg(S(:,4)));
hold on
plot(tspan./T, rad2deg(OM_filtered), 'k');

figure
om_filtered = movmean(S(:,5), [10000 10000]);
plot(tspan./T, rad2deg(S(:,5)));
hold on
plot(tspan./T, rad2deg(om_filtered), 'k');

figure
th_filtered = movmean(S(:,6), [10000 10000]);
plot(tspan./T, S(:,6));
hold on
plot(tspan./T, th_filtered, 'k');

%%
TH = TH';
