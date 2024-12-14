% LAB 5: chapeter 4
% 4 Dec 2024

% vectors are expressed in heliocentric ecliptic frame
% hypotesis: earth orbit is circular

clear all
close all
clc

mu_E = astroConstants(13);
mu_S = astroConstants(4);
AU = astroConstants(2);

v_inf_m = [15.1; 0; 0]; %[km/s]
imp_par = 9200; %[km]
r_E = AU*[1; 0; 0];
r_E_norm = norm(r_E);

% planetocentric hyperbola
v_inf = norm(v_inf_m); % v_inf = norm(v_inf_p) = norm(v_inf_p)
v_inf_m_dir = v_inf_m./v_inf; % direction of v_inf_min
a_HYP = -mu_E/v_inf^2;
delta = 2*atan2(-a_HYP, imp_par);
e_HYP = 1/sin(delta/2);
rP_HYP = a_HYP*(1-e_HYP);

% earth velocity
v_E_norm = sqrt(mu_S/r_E_norm);
v_E = [0; v_E_norm; 0]; % heliocentric - circular orbit

dV_norm = 2*v_inf*sin(delta/2);
beta = (pi-delta)/2;

v_pericenter = sqrt((mu_E*(1+e_HYP))/rP_HYP);

% leading
imp_par_vect_l = [0; imp_par; 0];
u_l = cross(imp_par_vect_l, v_inf_m)/norm(cross(imp_par_vect_l, v_inf_m));
v_inf_p_l = rotation_vector_Rodrigues(v_inf_m, u_l, delta);
V_m_l = v_E + v_inf_m;
V_p_l = v_E + v_inf_p_l;

tspan1 = [-20000000:1000: 0];
[T_ode, Y] = Orbit_Analysis(r_E, -V_m_l, mu_S, tspan1, 'non_perturbed', 'no_plot_orbit', 'no_plots');
tspan2 = [0:100: 20000000];
[T_ode2, Y2] = Orbit_Analysis(r_E, V_p_l, mu_S, tspan2, 'non_perturbed', 'no_plot_orbit', 'no_plots');


figure
plot3(Y(:,1), Y(:,2), Y(:, 3), 'b', LineWidth=2);
hold on
plot3(Y2(:,1), Y2(:,2), Y2(:, 3), 'r', LineWidth=2);
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
plot3(0, 0, 0, 'y*', LineWidth=3)
plot3(0, -1, 0, 'c*', LineWidth=3)
axis equal
grid on
legend('Before', 'After')
title('Heliocentric trajectories in HECI frame 3D')

% hyperbola 2d
a = -a_HYP; % Semiasse maggiore (valore esempio)
e = e_HYP; % Eccentricità (valore esempio)
c = a_HYP * e_HYP; % Distanza focale
b = a_HYP * sqrt(e_HYP^2 - 1); % Calcolo del semiasse minore

x = linspace(rP_HYP, r_E_norm/1000, 100000);
x_as_l = linspace(-r_E_norm/5000, c, 100000);

% equation of the hyp with center in (c2, 0)
y1 = b*sqrt((x+c).^2./a^2 - 1);
y2 = -b*sqrt((x+c).^2./a^2 - 1);

% asymptote
%as_l = @(X) -b/a*(X-c);

% 3D plot
Terra_3D
plot3(x, y2, zeros(size(x)), 'b', x, y1, zeros(size(x)), 'r', LineWidth=2);
hold on
plot3(0, 0, 0, 'c*', LineWidth=5);
%plot3(x_as_l, as_i(x_as_l), zeros(size(x_as_l)), 'm', linewidth=2);
%plot3(x_as_o, as_o(x_as_o), zeros(size(x_as_o)), 'g', linewidth=2);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title(' Flyby in Earth-centred frame parallel to HECI');
grid on;
axis equal;

%%
tspan = [0:1000: 100000000000];
[T_ode, Y_HYP] = Orbit_Analysis([-100000000; imp_par; 0], v_inf_m, mu_S, tspan, 'non_perturbed', 'no_plot_orbit', 'no_plots');

figure
plot3(Y_HYP(:,1), Y_HYP(:,2), Y_HYP(:, 3), 'b', LineWidth=2);
hold on
% plot3(0, AU, 0, 'c*', LineWidth=3)
xlabel('X [km]'); 
ylabel('Y [km]');
%zlabel('Z [km]');
axis equal

%%
% trailing
imp_par_vect_t = [0; -imp_par; 0];
u_t = cross(imp_par_vect_t, v_inf_m)/norm(cross(imp_par_vect_t, v_inf_m));
v_inf_p_t = rotation_vector_Rodrigues(v_inf_m, u_t, delta);
V_m_t = v_E + v_inf_m;
V_p_t = v_E + v_inf_p_t;

tspan1 = [-20000000:100: 0];
[T_ode3, Y3] = Orbit_Analysis(r_E, -V_m_t, mu_S, tspan1, 'non_perturbed', 'no_plot_orbit', 'no_plots');
tspan2 = [0:100: 20000000];
[T_ode4, Y4] = Orbit_Analysis(r_E, V_p_t, mu_S, tspan2, 'non_perturbed', 'no_plot_orbit', 'no_plots');

figure
plot3(Y3(:,1), Y3(:,2), Y3(:, 3), 'r', LineWidth=2);
hold on
plot3(Y4(:,1), Y4(:,2), Y4(:, 3), 'b', LineWidth=2);
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
plot3(0, 0, 0, 'y*', LineWidth=3)
axis equal


% under
imp_par_vect_u = [0; 0; -imp_par];
u_u = cross(imp_par_vect_u, v_inf_m)/norm(cross(imp_par_vect_u, v_inf_m));
v_inf_p_u = rotation_vector_Rodrigues(v_inf_m, u_u, delta);
V_m_u = v_E + v_inf_m;
V_p_u = v_E + v_inf_p_u;

tspan1 = [-20000000:100: 0];
[T_ode5, Y5] = Orbit_Analysis(r_E, -V_m_u, mu_S, tspan1, 'non_perturbed', 'no_plot_orbit', 'no_plots');
tspan2 = [0:100: 20000000];
[T_ode6, Y6] = Orbit_Analysis(r_E, V_p_u, mu_S, tspan2, 'non_perturbed', 'no_plot_orbit', 'no_plots');

figure
plot3(Y5(:,1), Y5(:,2), Y5(:, 3), 'r', LineWidth=2);
hold on
plot3(Y6(:,1), Y6(:,2), Y6(:, 3), 'b', LineWidth=2);
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
plot3(0, 0, 0, 'y*', LineWidth=3)
axis equal


%% 

mu_E = astroConstants(13);
mu_S = astroConstants(4);
AU = astroConstants(2);

v_inf_m = [15.1; 0; 0]; %[km/s]
imp_par = linspace(9200, 13200, 5); %[km]
r_E = AU*[1; 0; 0];
r_E_norm = norm(r_E);

% planetocentric hyperbola
v_inf = norm(v_inf_m); % v_inf = norm(v_inf_p) = norm(v_inf_p)
v_inf_m_dir = v_inf_m./v_inf; % direction of v_inf_min
a_HYP = -mu_E/v_inf^2;
delta = 2*atan2(-a_HYP, imp_par);
e_HYP = 1./sin(delta/2);
rP_HYP = a_HYP.*(1-e_HYP);

% earth velocity
v_E_norm = sqrt(mu_S/r_E_norm);
v_E = [0; v_E_norm; 0]; % heliocentric - circular orbit

dV_norm = 2*v_inf*sin(delta./2);
beta = (pi-delta)/2;

%v_pericenter = sqrt((mu_E*(1+e_HYP))/rP_HYP);

% leading
figure
u_l = -cross(r_E, v_E)/norm(cross(r_E, v_E));
v_inf_p_l = rotation_vector_Rodrigues(v_inf_m, u_l, delta);
V_m_l = v_E + v_inf_m;

tspan1 = [-20000000:100: 0];
[T_ode, Y] = Orbit_Analysis(r_E, -V_m_l, mu_S, tspan1, 'non_perturbed', 'no_plot_orbit', 'no_plots');


plot3(Y(:,1), Y(:,2), Y(:, 3), 'r', LineWidth=2);
hold on
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
plot3(0, 0, 0, 'y*', LineWidth=3)
axis equal
grid on

for i = 1:5
    V_p_l = v_E + v_inf_p_l(:,i);
    %V_p_l = V_p_l_v(:,i);
    tspan2 = [0:100: 20000000];
    [T_ode2, Y2] = Orbit_Analysis(r_E, V_p_l, mu_S, tspan2, 'non_perturbed', 'no_plot_orbit', 'no_plots');
    plot3(Y2(:,1), Y2(:,2), Y2(:, 3), 'b', LineWidth=2);
end

%% Powered gravity assist
clear all 
close all
clc

%

% CONSTANTS
mu_E = astroConstants(13);
mu_S = astroConstants(4);
AU = astroConstants(2);

% DATA
V_m = [31.5; 5.2; 0.0]; %[km/s]
V_p = [36.0; 0.0; 0.0]; %[km/2]
r_E = AU*[0; -1; 0];

% solve to get the v_inf + and - and delta
r_E_norm = norm(r_E);

v_E_dir = cross([0; 0; 1], [0; -1; 0]);
v_E_norm = sqrt(mu_S/r_E_norm);
v_E = v_E_norm*v_E_dir;

v_inf_p = V_p - v_E;
v_inf_m = V_m - v_E;

v_inf_m_norm = norm(v_inf_m);
v_inf_p_norm = norm(v_inf_p);

a_m = -mu_E/v_inf_m_norm^2;
a_p = -mu_E/v_inf_p_norm^2;

delta = acos(dot(v_inf_m, v_inf_p)/(v_inf_m_norm*v_inf_p_norm));
delta_deg = rad2deg(delta);

% solve the equation: delta = f(rp, v_inf_p_norm, v_inf_p_norm) = delta_m/2 + delta_p/2
em = @(rP) 1 + rP*v_inf_m_norm^2/mu_E;
ep = @(rP) 1 + rP*v_inf_p_norm^2/mu_E;

delta_m = @(rP) 2*asin(1/em(rP));
delta_p = @(rP) 2*asin(1/ep(rP));

fun = @(rP) delta_m(rP)/2 + delta_p(rP)/2 -delta;

r_P_norm = fzero(fun, 50000); % sensitive for the starting condition but fsolve is even worst
                         % tried with different initial condition, you need
                         % to be close enough, expeccialy for minor initial
                         % conditions: [4346, up to 50000 is still fine and
                         % gives the sam results probaly even higher (but
                         % lower than 100000]

% eccentricities
e_m = em(r_P_norm);
e_p = ep(r_P_norm);

% the following velocities are all parallel, therefore you just have to
% work with their norms
% calculate the deltaV at pericenter
v_P_m_norm = sqrt((mu_E*(1+e_m))/r_P_norm);
v_P_p_norm = sqrt((mu_E*(1+e_p))/r_P_norm);
dv_P = abs(v_P_p_norm-v_P_m_norm);

% deltaV fly-by = norm of the difference of the v_inf
dV_fb = norm(v_inf_p - v_inf_m);

% PLOTS 
% Heliocentric trajectories in HECI frame
tspan1 = [-20000000:1000: 0];
[T_ode, Y] = Orbit_Analysis(r_E, -V_m, mu_S, tspan1, 'non_perturbed', 'no_plot_orbit', 'no_plots');
tspan2 = [0:100: 20000000];
[T_ode2, Y2] = Orbit_Analysis(r_E, V_p, mu_S, tspan2, 'non_perturbed', 'no_plot_orbit', 'no_plots');

% 3D
figure
plot3(Y(:,1)./AU, Y(:,2)./AU, Y(:, 3)./AU, 'b', LineWidth=2);
hold on
plot3(Y2(:,1)./AU, Y2(:,2)./AU, Y2(:, 3)./AU, 'r', LineWidth=2);
xlabel('X [AU]'); 
ylabel('Y [AU]');
zlabel('Z [AU]');
plot3(0, 0, 0, 'y*', LineWidth=10)
plot3(0, -1, 0, 'c*', LineWidth=3)
axis equal
grid on
legend('Before', 'After')
title('Heliocentric trajectories in HECI frame 3D')

% 2D
figure
plot(Y(:,1)./AU, Y(:,2)./AU, 'b', LineWidth=2);
hold on
plot(Y2(:,1)./AU, Y2(:,2)./AU, 'r', LineWidth=2);
xlabel('X [AU]'); 
ylabel('Y [AU]');
plot(0, 0, 'y*', LineWidth=10)
plot(0, -1, 'c*', LineWidth=3)
axis equal
grid on
legend('Before', 'After')
title('Heliocentric trajectories in HECI frame 2D')


% incoming hyperbola
a = -a_m; % Semiasse maggiore (valore esempio)
e = e_m; % Eccentricità (valore esempio)
c = a * e; % Distanza focale
b = a * sqrt(e^2 - 1); % Calcolo del semiasse minore

x = linspace(-r_E_norm/5000, r_P_norm, 100000);
x_as_i = linspace(-r_E_norm/5000, c, 100000);

% equation of the hyp with center in (c, 0)
y1_i = b*sqrt((x-c).^2./a^2 - 1);
y2_i = -b*sqrt((x-c).^2./a^2 - 1);

% asymptote
as_i = @(X) b/a*(X-c);


% outcoming hyperbola
a_o = -a_p; % Semiasse maggiore (valore esempio)
e_o = e_p; % Eccentricità (valore esempio)
c_o = a_o * e_o; % Distanza focale
b_o = a_o * sqrt(e_o^2 - 1); % Calcolo del semiasse minore

x_as_o = linspace(-r_E_norm/5000, c_o, 100000);

% equation of the hyp with center in (c2, 0)
y1_o = b_o*sqrt((x-c_o).^2./a_o^2 - 1);
y2_o = -b_o*sqrt((x-c_o).^2./a_o^2 - 1);

% asymptote
as_o = @(X) -b_o/a_o*(X-c_o);

% plot
% plot
figure
plot(x, y2_i, 'b', LineWidth=2);
hold on
plot(x_as_i, as_i(x_as_i), 'm', linewidth=2);
plot(0, 0, 'c*', LineWidth=5);
grid on
axis equal
plot(x, y1_o, 'r', LineWidth=2);
plot(x_as_o, as_o(x_as_o), 'g', linewidth=2);
plot(0, 0, 'c*', LineWidth=5);
title('Plot in flyby perifocal frame')

% 3D plot - ORIENTAZIONE NEL PIANO NON VERITIERA
Terra_3D
plot3(x, y2_i, zeros(size(x)), 'b', LineWidth=2);
hold on
plot3(x, y1_o, zeros(size(x)), 'r', LineWidth=2);
plot3(0, 0, 0, 'c*', LineWidth=5);
plot3(x_as_i, as_i(x_as_i), zeros(size(x_as_i)), 'm', linewidth=2);
plot3(x_as_o, as_o(x_as_o), zeros(size(x_as_o)), 'g', linewidth=2);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title(' Flyby in Earth-centred frame parallel to HECI');
grid on;
axis equal;

%
v_inf_m_dir = 1/v_inf_m_norm * v_inf_m;
v_inf_p_dir = 1/v_inf_p_norm * v_inf_p;
beta_m = pi/2 - delta_m(r_P_norm)/2;
beta_p = pi/2 - delta_p(r_P_norm)/2;
apse_line_dir = -rotation_vector_Rodrigues(v_inf_p_dir, [0; 0; 1],-beta_p); % you want the oppositie of it to decompose c
cx = c*apse_line_dir(1);
cy = c*apse_line_dir(2);
c_ox = c_o*apse_line_dir(1);
c_oy = c_o*apse_line_dir(2);

% asymptotes and apse line
as_m = @(X) v_inf_m_dir(2)/v_inf_m_dir(1) * (X-cx) + cy;
as_p = @(X) v_inf_p_dir(2)/v_inf_p_dir(1) * (X-c_ox) + c_oy;
apse_line = @(X) apse_line_dir(2)/apse_line_dir(1) * X; 

% rotation 
rot_angle = acos(dot([1,0,0], apse_line_dir)/(norm(apse_line_dir)));
% Applicazione della matrice di rotazione
R = [cos(rot_angle), -sin(rot_angle); sin(rot_angle), cos(rot_angle)];

% incoming hyperbola after rotation
rotated_y1_i = R * [x; y1_i];
rotated_y2_i = R * [x; y2_i];

x_rotated1_i = rotated_y1_i(1, :);
y_rotated1_i = rotated_y1_i(2, :);
x_rotated2_i = rotated_y2_i(1, :);
y_rotated2_i = rotated_y2_i(2, :);


% upcoming hyperbola after rotation
rotated_y1_o = R * [x; y1_o];
rotated_y2_o = R * [x; y2_o];

x_rotated1_o = rotated_y1_o(1, :);
y_rotated1_o = rotated_y1_o(2, :);
x_rotated2_o = rotated_y2_o(1, :);
y_rotated2_o = rotated_y2_o(2, :);

% 3D plot with the correct orientation of the asymptote and the hyperbola
Terra_3D
hold on
plot3([-r_E_norm/2000:1:r_E_norm/2000], as_m([-r_E_norm/2000:1:r_E_norm/2000]),0*[-r_E_norm/2000:1:r_E_norm/2000], 'm');
plot3([-r_E_norm/5000:1:r_E_norm/5000], as_p([-r_E_norm/5000:1:r_E_norm/5000]),0*[-r_E_norm/5000:1:r_E_norm/5000], 'g');
plot3([-r_E_norm/5000:1:r_E_norm/5000], apse_line([-r_E_norm/5000:1:r_E_norm/5000]),0*[-r_E_norm/5000:1:r_E_norm/5000], 'k');
plot3(x_rotated1_i, y_rotated1_i, 0*x_rotated1_i, 'b') %, x_rotated2_i, y_rotated2_i, 'b');
hold on
plot3(x_rotated2_o, y_rotated2_o, 0*x_rotated1_i, 'r') %, x_rotated2_i, y_rotated2_i, 'b');
axis equal
grid on
title('Flyby in Earth-centred frame parallel to HECI');

%%
% 3D plot - ORIENTAZIONE NEL PIANO NON VERITIERA
Terra_3D
hold on
plot3(x_as_i, as_i(x_as_i), zeros(size(x_as_i)), 'm', linewidth=2);
plot3(x_as_o, as_o(x_as_o), zeros(size(x_as_o)), 'g', linewidth=2);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title(' Flyby in Earth-centred frame parallel to HECI');
grid on;
axis equal;