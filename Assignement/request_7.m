%% 7TH REQUEST
% two objects (non operative satellites) orbiting around Earth
% Data from Space Track (TLE) and NASA Horizon System (ephemerides)

clear all
close all
clc

% Data for Earth as primary planet
mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 1 1 0 0 0]);

%% IMAGE (High Elliptical Orbit) 
% upload the information from IMAGEeph by NASA/JPL Horizon
load('IMAGEeph.mat');

a_eph_image = IMAGEeph(:,10);
e_eph_image = IMAGEeph(:,1);
i_eph_image = deg2rad(IMAGEeph(:,3));
OM_eph_image = deg2rad(IMAGEeph(:,4));
om_eph_image = deg2rad(IMAGEeph(:,5));
th_eph_image = deg2rad(IMAGEeph(:,9));

T_orb_eph_image = 2*pi*sqrt(a_eph_image(1)^3/mu_E);

for i = 1:length(a_eph_image)
    [rr, vv] = par2car(a_eph_image(i), e_eph_image(i), i_eph_image(i), OM_eph_image(i), om_eph_image(i), th_eph_image(i), mu_E);
    r_image(i, :) = rr;
    v_image(i, :) = vv;
end

Terra_3D
hold on
plot3(r_image(:,1), r_image(:,2), r_image(:,3))
grid on 
axis equal
title('ephemerides plot')

% propagation of the orbit
tspan3 = [0:60*4:60*24*60*60];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
eph_0_image = [a_eph_image(1); e_eph_image(1); i_eph_image(1); OM_eph_image(1); om_eph_image(1); th_eph_image(1)];

tic
[T_pert_eph_image, Y_pert_eph_image] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_image, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph_image, 1)
    [rr, vv] = par2car(Y_pert_eph_image(i,1), Y_pert_eph_image(i,2), Y_pert_eph_image(i,3), Y_pert_eph_image(i,4), Y_pert_eph_image(i,5), Y_pert_eph_image(i,6), mu_E);
    rrr_image(i, :) = rr;
    vvv_image(i, :) = vv;
end

Terra_3D
hold on
scatter3(rrr_image(:,1), rrr_image(:,2), rrr_image(:,3), 3, T_pert_eph_image, 'filled');
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';

% Behaviour
% unwrap of the angles
OM_eph_image = unwrap(OM_eph_image);
om_eph_image = unwrap(om_eph_image);
th_eph_image = unwrap(th_eph_image);

% to have the behaviour expressed in days
vect = tspan3./(60*60*24);

% plot of the behavior of the keplerian elements to compare
figure
subplot(3,2,1)
plot(vect, a_eph_image(:));
hold on
plot(vect, Y_pert_eph_image(:,1));
grid on
title('Semi-major axis evolution over time', interpreter='latex', FontSize=13, FontWeight='bold')
legend('ephemerides', 'propagator', interpreter='latex', FontSize=12)
ylabel('a [km]', interpreter='latex', FontSize=12)
xlabel('Time [day]', interpreter='latex', FontSize=1)

subplot(3,2,2)
plot(vect, e_eph_image(:));
hold on
plot(vect, Y_pert_eph_image(:,2));
grid on
title('Eccentricity evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('e [-]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,3)
plot(vect, rad2deg(i_eph_image(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_image(:,3)));
grid on
title('Inclination evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('i [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,4)
plot(vect, rad2deg(OM_eph_image(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_image(:,4)));
grid on
title('RAAN evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\Omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,5)
plot(vect, rad2deg(om_eph_image(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_image(:,5)));
grid on
title('Argument of pericenter evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')


subplot(3,2,6)
plot(vect,  rad2deg(th_eph_image(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_image(:,6)));
grid on
title('True anomaly evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
xlabel('Time [day]', interpreter='latex')
ylabel('$\theta$ [deg]', interpreter='latex')

%% --- Low Earth Orbit
% upload the information from DMSP13_eph
clear all
close all
clc

% Data for Earth as primary planet
mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 12 14 0 0 0]);
%%
load('DMSPeph.mat');
%%
a_eph_DMSP = DMSPeph(:,10);
e_eph_DMSP = DMSPeph(:,1);
i_eph_DMSP = deg2rad(DMSPeph(:,3));
OM_eph_DMSP = deg2rad(DMSPeph(:,4));
om_eph_DMSP = deg2rad(DMSPeph(:,5));
th_eph_DMSP = deg2rad(DMSPeph(:,9));

T_orb_eph_DMSP = 2*pi*sqrt(a_eph_DMSP(1)^3/mu_E);

tspan3 = [0:60*4:31*24*60*60];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

eph_0_DMSP = [a_eph_DMSP(1); e_eph_DMSP(1); i_eph_DMSP(1); OM_eph_DMSP(1); om_eph_DMSP(1); th_eph_DMSP(1)];

for i = 1:length(a_eph_DMSP)
    [rr, vv] = par2car(a_eph_DMSP(i), e_eph_DMSP(i), i_eph_DMSP(i), OM_eph_DMSP(i), om_eph_DMSP(i), th_eph_DMSP(i), mu_E);
    r_DMSP(i, :) = rr;
    v_DMSP(i, :) = vv;
end

Terra_3D
hold on
plot3(r_DMSP(:,1), r_DMSP(:,2), r_DMSP(:,3))
grid on 
axis equal

tic
[T_pert_eph_DMSP, Y_pert_eph_DMSP] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_DMSP, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph_DMSP, 1)
    [rr, vv] = par2car(Y_pert_eph_DMSP(i,1), Y_pert_eph_DMSP(i,2), Y_pert_eph_DMSP(i,3), Y_pert_eph_DMSP(i,4), Y_pert_eph_DMSP(i,5), Y_pert_eph_DMSP(i,6), mu_E);
    rrr_DMSP(i, :) = rr;
    vvv_DMSP(i, :) = vv;
end

Terra_3D
hold on
scatter3(rrr_DMSP(:,1), rrr_DMSP(:,2), rrr_DMSP(:,3), 3, T_pert_eph_DMSP, 'filled');
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';

%% Behaviour

OM_eph_DMSP = unwrap(OM_eph_DMSP);
om_eph_DMSP = unwrap(om_eph_DMSP);
th_eph_DMSP = unwrap(th_eph_DMSP);

vect = tspan3./(60*60*24);
figure
subplot(3,2,1)
plot(vect, a_eph_DMSP(:));
hold on
plot(vect, Y_pert_eph_DMSP(:,1));
grid on
xlim([0 31])
ylim([7180 7260])
title('Semi-major axis evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('a [km]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,2)
plot(vect, e_eph_DMSP(:));
hold on
plot(vect, Y_pert_eph_DMSP(:,2));
grid on
xlim([0 31])
ylim([-2e-3 6e-3])
title('Eccentricity evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('e [-]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,3)
plot(vect, rad2deg(i_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,3)));
grid on
xlim([0 31])
ylim([98.7 99.2])
title('Inclination evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('i [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,4)
plot(vect, rad2deg(OM_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,4)));
grid on
xlim([0 31])
title('RAAN evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\Omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,5)
plot(vect, rad2deg(om_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,5)));
grid on
xlim([0 31])
ylim([-50 300])
title('Argument of pericenter evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')


subplot(3,2,6)
plot(vect,  rad2deg(th_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,6)));
grid on
xlim([0 30])
title('True anomaly evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
xlabel('Time [day]', interpreter='latex')
ylabel('$\theta$ [deg]', interpreter='latex')

%%

%% --- Low Earth Orbit
% upload the information from DMSP13_eph
clear all
close all
clc

% Data for Earth as primary planet
mu_E = astroConstants(13);
J2_E = astroConstants(33);
R_E = astroConstants(23);
eps_E = astroConstants(63);
start_date_2 = date2mjd2000([2024 12 14 0 0 0]);
%%
load('DMSPeph.mat');
%%
a_eph_DMSP = DMSPeph(:,10);
e_eph_DMSP = DMSPeph(:,1);
i_eph_DMSP = deg2rad(DMSPeph(:,3));
OM_eph_DMSP = deg2rad(DMSPeph(:,4));
om_eph_DMSP = deg2rad(DMSPeph(:,5));
th_eph_DMSP = deg2rad(DMSPeph(:,9));

T_orb_eph_DMSP = 2*pi*sqrt(a_eph_DMSP(1)^3/mu_E);

tspan3 = [0:60*4:31*24*60*60];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

eph_0_DMSP = [a_eph_DMSP(1); e_eph_DMSP(1); i_eph_DMSP(1); OM_eph_DMSP(1); om_eph_DMSP(1); th_eph_DMSP(1)];

for i = 1:length(a_eph_DMSP)
    [rr, vv] = par2car(a_eph_DMSP(i), e_eph_DMSP(i), i_eph_DMSP(i), OM_eph_DMSP(i), om_eph_DMSP(i), th_eph_DMSP(i), mu_E);
    r_DMSP(i, :) = rr;
    v_DMSP(i, :) = vv;
end

Terra_3D
hold on
plot3(r_DMSP(:,1), r_DMSP(:,2), r_DMSP(:,3))
grid on 
axis equal

tic
[T_pert_eph_DMSP, Y_pert_eph_DMSP] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_DMSP, options);
elapsedtime_gauss = toc;
fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)

for i = 1:size(Y_pert_eph_DMSP, 1)
    [rr, vv] = par2car(Y_pert_eph_DMSP(i,1), Y_pert_eph_DMSP(i,2), Y_pert_eph_DMSP(i,3), Y_pert_eph_DMSP(i,4), Y_pert_eph_DMSP(i,5), Y_pert_eph_DMSP(i,6), mu_E);
    rrr_DMSP(i, :) = rr;
    vvv_DMSP(i, :) = vv;
end

Terra_3D
hold on
scatter3(rrr_DMSP(:,1), rrr_DMSP(:,2), rrr_DMSP(:,3), 3, T_pert_eph_DMSP, 'filled');
colormap('default')
colorbar
c = colorbar;
c.Label.String = 'Time [s]';

%% Behaviour

OM_eph_DMSP = unwrap(OM_eph_DMSP);
om_eph_DMSP = unwrap(om_eph_DMSP);
th_eph_DMSP = unwrap(th_eph_DMSP);

vect = tspan3./(60*60*24);
figure
subplot(3,2,1)
plot(vect, a_eph_DMSP(:));
hold on
plot(vect, Y_pert_eph_DMSP(:,1));
grid on
xlim([0 31])
ylim([7180 7260])
title('Semi-major axis evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('a [km]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,2)
plot(vect, e_eph_DMSP(:));
hold on
plot(vect, Y_pert_eph_DMSP(:,2));
grid on
xlim([0 31])
ylim([-2e-3 6e-3])
title('Eccentricity evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('e [-]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,3)
plot(vect, rad2deg(i_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,3)));
grid on
xlim([0 31])
ylim([98.7 99.2])
title('Inclination evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('i [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,4)
plot(vect, rad2deg(OM_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,4)));
grid on
xlim([0 31])
title('RAAN evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\Omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')

subplot(3,2,5)
plot(vect, rad2deg(om_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,5)));
grid on
xlim([0 31])
ylim([-50 300])
title('Argument of pericenter evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
ylabel('$\omega$ [deg]', interpreter='latex')
xlabel('Time [day]', interpreter='latex')


subplot(3,2,6)
plot(vect,  rad2deg(th_eph_DMSP(:)));
hold on
plot(vect, rad2deg(Y_pert_eph_DMSP(:,6)));
grid on
xlim([0 30])
title('True anomaly evolution over time', interpreter='latex')
legend('ephemerides', 'propagator', interpreter='latex')
xlabel('Time [day]', interpreter='latex')
ylabel('$\theta$ [deg]', interpreter='latex')
%%
% %% Envisat (Low Earth Orbit) 
% % upload the information from envisat_eph
% %load('envisateph.mat');
% load('envisateph2.mat');
% 
% % a_eph_envisat = envisateph(:,10);
% % e_eph_envisat = envisateph(:,1);
% % i_eph_envisat = deg2rad(envisateph(:,3));
% % OM_eph_envisat = deg2rad(envisateph(:,4));
% % om_eph_envisat = deg2rad(envisateph(:,5));
% % th_eph_envisat = deg2rad(envisateph(:,9));
% 
% a_eph_envisat = envisateph2(:,10);
% e_eph_envisat = envisateph2(:,1);
% i_eph_envisat = deg2rad(envisateph2(:,3));
% OM_eph_envisat = deg2rad(envisateph2(:,4));
% om_eph_envisat = deg2rad(envisateph2(:,5));
% th_eph_envisat = deg2rad(envisateph2(:,9));
% 
% T_orb_eph_envisat = 2*pi*sqrt(a_eph_envisat(1)^3/mu_E);
% nn = 864000/T_orb_eph_envisat;
% tspan3 = [0:60*2:31*24*60*60];
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% 
% eph_0_envisat = [a_eph_envisat(1); e_eph_envisat(1); i_eph_envisat(1); OM_eph_envisat(1); om_eph_envisat(1); th_eph_envisat(1)];
% 
% for i = 1:length(a_eph_envisat)
%     [rr, vv] = par2car(a_eph_envisat(i), e_eph_envisat(i), i_eph_envisat(i), OM_eph_envisat(i), om_eph_envisat(i), th_eph_envisat(i), mu_E);
%     r_envisat(i, :) = rr;
%     v_envisat(i, :) = vv;
% end
% 
% Terra_3D
% hold on
% plot3(r_envisat(:,1), r_envisat(:,2), r_envisat(:,3))
% grid on 
% axis equal
% 
% tic
% [T_pert_eph_envisat, Y_pert_eph_envisat] = ode113(@(t, k_el) eq_motion_Gauss_RSW(t, k_el, @(t, k_el) acc_pert_function(t, k_el, J2_E, mu_E, R_E, start_date_2, eps_E, 3), mu_E), tspan3, eph_0_envisat, options);
% elapsedtime_gauss = toc;
% fprintf ('Gaussian simulation from eph: %.4f s. \n', elapsedtime_gauss)
% 
% for i = 1:size(Y_pert_eph_envisat, 1)
%     [rr, vv] = par2car(Y_pert_eph_envisat(i,1), Y_pert_eph_envisat(i,2), Y_pert_eph_envisat(i,3), Y_pert_eph_envisat(i,4), Y_pert_eph_envisat(i,5), Y_pert_eph_envisat(i,6), mu_E);
%     rrr_envisat(i, :) = rr;
%     vvv_envisat(i, :) = vv;
% end
% 
% Terra_3D
% hold on
% scatter3(rrr_envisat(:,1), rrr_envisat(:,2), rrr_envisat(:,3), 3, T_pert_eph_envisat, 'filled');
% colormap('default')
% colorbar
% c = colorbar;
% c.Label.String = 'Time [s]';
% 
% %% Behaviour
% 
% OM_eph_envisat = unwrap(OM_eph_envisat);
% om_eph_envisat = unwrap(om_eph_envisat);
% th_eph_envisat = unwrap(th_eph_envisat);
% 
% vect = tspan3./(60*60*24);
% figure
% subplot(3,2,1)
% plot(vect, a_eph_envisat(:));
% hold on
% plot(vect, Y_pert_eph_envisat(:,1));
% grid on
% title('Semi-major axis evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% ylabel('a [km]', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% 
% subplot(3,2,2)
% plot(vect, e_eph_envisat(:));
% hold on
% plot(vect, Y_pert_eph_envisat(:,2));
% grid on
% title('Eccentricity evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% ylabel('e [-]', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% 
% subplot(3,2,3)
% plot(vect, rad2deg(i_eph_envisat(:)));
% hold on
% plot(vect, rad2deg(Y_pert_eph_envisat(:,3)));
% grid on
% title('Inclination evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% ylabel('i [deg]', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% 
% subplot(3,2,4)
% plot(vect, rad2deg(OM_eph_envisat(:)));
% hold on
% plot(vect, rad2deg(Y_pert_eph_envisat(:,4)));
% grid on
% title('RAAN evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% ylabel('$\Omega$ [deg]', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% 
% subplot(3,2,5)
% plot(vect, rad2deg(om_eph_envisat(:)));
% hold on
% plot(vect, rad2deg(Y_pert_eph_envisat(:,5)));
% grid on
% title('Argument of pericenter evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% ylabel('$\omega$ [deg]', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% 
% 
% subplot(3,2,6)
% plot(vect,  rad2deg(th_eph_envisat(:)));
% hold on
% plot(vect, rad2deg(Y_pert_eph_envisat(:,6)));
% grid on
% title('True anomaly evolution over time', interpreter='latex')
% legend('ephemerides', 'propagator', interpreter='latex')
% xlabel('Time [day]', interpreter='latex')
% ylabel('$\theta$ [deg]', interpreter='latex')

