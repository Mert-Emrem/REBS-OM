%% LAB 4 
% 13/11/2024

%% ES 1 
clear all
close all
clc

mu = astroConstants(13);
r1 = [-21800; 37900; 0];   %km
r2 = [27300; 27700; 0];    %km
dt = 15*60*60 + 6*60 + 40; %s
[a_T, p_T, e_T, error_flag, v1, v2, t_par, theta] = lambertMR(r1, r2, dt, mu, 0, 0, 0);
v1=v1';
T = 2*pi*sqrt(a_T^3/mu);
tspan = [1:0.1:T];
[T_ode, Y] = Orbit_Analysis(r1, v1, mu, tspan, 'non_perturbed', 'plot_orbit', 'no_plots');

hold on 
plot3(r1(1), r1(2), r1(3), 'r*', LineWidth=4)
plot3(r2(1), r2(2), r2(3), 'r*', LineWidth=4)
quiver3(0, 0, 0, r1(1), r1(2), r1(3), 0);
quiver3(0, 0, 0, r2(1), r2(2), r2(3), 0);

%% E2 2
clear all
close all
clc

mu = astroConstants(13);

a1 = 12500;
e1 = 0;
i1 = 0;
OMEGA1 = 0;
omega1 = 0;
theta1 = deg2rad(120);

a2 = 9500;
e2 = 0.3;
i2 = 0;
OMEGA2 = 0;
omega2 = 0; 
theta2 = deg2rad(250);

tof = 3300;
[r1, v1] = par2car(a1, e1, i1, OMEGA1, omega1, theta1, mu);
[r2, v2] = par2car(a2, e2, i2, OMEGA2, omega2, theta2, mu);

[a_T, p_T, e_T, error_flag, v1_T, v2_T, t_par, theta] = lambertMR(r1, r2, tof, mu, 0, 0, 0);
v1_T = v1_T';
v2_T = v2_T';
dV1 = norm(v1_T - v1);
dV2 = norm(v2 - v2_T);
dVtot = dV1 + dV2;

tspan = [0:0.1:tof];
t1 = [0:0.1:2*pi*sqrt(a1^3/mu)];
t2 = [0:0.1:2*pi*sqrt(a2^3/mu)];
[~, YT] = Orbit_Analysis(r1, v1_T, mu, tspan, 'non_perturbed', 'no_plot_orbit', 'no_plots');
[~, Y1] = Orbit_Analysis(r1, v1, mu, t1, 'non_perturbed', 'no_plot_orbit', 'no_plots');
[~, Y2] = Orbit_Analysis(r2, v2, mu, t2, 'non_perturbed', 'no_plot_orbit', 'no_plots');

Terra_3D
hold on
plot3(Y1(:,1), Y1(:,2), Y1(:, 3), 'r', LineWidth=2);
plot3(Y2(:,1), Y2(:,2), Y2(:, 3), 'b', LineWidth=2);
plot3(YT(:,1), YT(:,2), YT(:, 3), 'g', LineWidth=2);
plot3(r1(1), r1(2), r1(3), 'k*', LineWidth=4)
plot3(r2(1), r2(2), r2(3), 'k*', LineWidth=4)
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');

%% es 3
clear all
close all
clc

time1 = date2mjd2000([2003, 4, 1, 0, 0, 0]);
time2 = date2mjd2000([2003, 8, 1, 23, 59, 59]);

time3 = date2mjd2000([2003, 9, 1, 0, 0, 0]);
time4 = date2mjd2000([2004, 3, 1, 23, 59, 59]);

fracDay = hms2fracday(12, 0, 0);
window_launch = [time1:fracDay:time2];
window_arrival = [time3:fracDay:time4];

TIME = zeros(length(window_launch), length(window_arrival));
DV = zeros(length(window_launch), length(window_arrival));
VL = zeros(length(window_launch), 1);
VA = zeros(length(window_arrival), 1);
V1T = zeros(length(window_launch), length(window_arrival));
V2T = zeros(length(window_launch), length(window_arrival));

for i = 1:1:length(window_launch)
    L = window_launch(i);
    [kepEl,ksun] = uplanet(L, 3);
    [rEl, vEl] = par2car(kepEl(1), kepEl(2), kepEl(3), kepEl(4), kepEl(5), kepEl(6), ksun); 
    l = L*24*3600;
    for j = 1:1:length(window_arrival)
        A = window_arrival(j);
        a = A*24*3600;
        dt = a-l;
        TIME(i,j) = dt;
        [kepMa,~] = uplanet(A, 4);
        [rMa, vMa] = par2car(kepMa(1), kepMa(2), kepMa(3), kepMa(4), kepMa(5), kepMa(6), ksun);  
        [a_T, ~, e_T, ~, v1_T, v2_T, t_par, theta] = lambertMR(rEl, rMa, dt, ksun, 0, 0, 0);
        v1_T = v1_T';
        v2_T = v2_T';
        dV1 = norm(v1_T - vEl);
        dV2 = norm(vMa - v2_T);
        dVtot = dV1 + dV2;
        DV(i, j) = dVtot;
    end
end


fracDay_xfun = hms2fracday(0, 0, 1);
window_launch_xfun = [time1:fracDay:time2];
window_arrival_xfun = [time3:fracDay:time4];
%DV_fun = interp2(window_launch', window_arrival', DV',  window_arrival_xfun', window_launch_xfun);
%[dv_opt] = min(DV, [], 'all');
%M_L = fminunc(DV_fun, dv_opt);

[dv_opt] = min(DV, [], 'all');
[ii, jj] = find (DV == dv_opt);
date_launch = window_launch(ii);
date_arrival = window_arrival(jj);
dt_transfer = TIME(ii, jj);
LAUNCH = mjd20002date(date_launch);
ARRIVAL = mjd20002date(date_arrival);

levels = [6, 7, 8, 9, 10];
figure
contour(window_launch, window_arrival, DV', levels);
colorbar
grid on
title('Porkchop plot');
xlabel('Launch');
ylabel('Arrival');


tspanE = [0:1000:365*24*3600];
tspanM = [0:1000:687*24*3600];
tspanT = [0:1000:dt_transfer];
[kepEl,ksun] = uplanet(date_launch, 3);
[rEL, vEL] = par2car(kepEl(1), kepEl(2), kepEl(3), kepEl(4), kepEl(5), kepEl(6), ksun);
[~, YE] = Orbit_Analysis(rEL, vEL, ksun, tspanE, 'non_perturbed', 'no_plot_orbit', 'no_plots');
[kepMa,~] = uplanet(A, 4);
[rMA, vMA] = par2car(kepMa(1), kepMa(2), kepMa(3), kepMa(4), kepMa(5), kepMa(6), ksun);  
[~, YM] = Orbit_Analysis(rMA, vMA, ksun, tspanM, 'non_perturbed', 'no_plot_orbit', 'no_plots');
[a_T, ~, e_T, ~, v1_T, v2_T, t_par, theta] = lambertMR(rEL, rMA, dt_transfer, ksun, 0, 0, 0);   
v1_T = v1_T';
[~, YT] = Orbit_Analysis(rEL, v1_T, ksun, tspanT, 'non_perturbed', 'no_plot_orbit', 'no_plots');

% Crea un vettore di date


% Crea alcuni dati casuali da plottare
figure
plot3(YE(:,1), YE(:,2), YE(:,3), 'b', LineWidth=2);
hold on
plot3(YM(:,1), YM(:,2), YM(:,3), 'r', LineWidth=2);
plot3(YT(:,1), YT(:,2), YT(:,3), 'g', LineWidth=2);
grid on



% Aggiungi griglia



