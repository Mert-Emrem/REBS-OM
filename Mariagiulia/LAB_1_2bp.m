%% LAB 1 
% done by myself on the 12th Oct 24 
% reference slide: "Lab Chapter 1: The two body problem"

%% EX1A_EARTH_2_BODY_PB

clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [26578.137; 0; 0]; %[km]
v0 = [0; 2.221; 3.173]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
%year = 365*24*60*60;
tspan = linspace(0, 2*T, 1000); % span for the ODE resolution

% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); % call this when you use the ode to set your own tollerances (adequate for the precision required for orbital porbelms

% perform the integration
[T_ode, Y] = ode113(@(t,y) ode_2pb(t, y, mu_E), tspan, y0, options);
                
                % choice of 113: more efficient than 45 with stingent tolleances
                % put the variables and call the system of equations
                % give the interval of time: t0 ti in which you want the evaluation tf
                % give the state vector 0
                % call the options about the tollerances
                % T = time, Y = matrix of the coordinatez x y z of the position at each time step asked

% Plot the results of the orbit
figure (1)
plot3(Y(:,1), Y(:,2), Y(:, 3), '-');
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit - Earth');
axis equal;
grid on;

% check the constants
r = Y(:, 1:3)';
v = Y(:, 4:6)';

% angular momentum
h_vect_0 = cross(r0, v0);
% h_vect_t = cross(Y(:,1:3)',Y(:,4:6)'); 
h_vect = cross(r, v);
h_norm = vecnorm(h_vect);

% energy
% en = dot(Y(:,4:6)', Y(:,4:6)')./2 - mu_E./vecnorm(Y(:,1:3)');
en = vecnorm(v).^2./2 - mu_E./vecnorm(r);

% eccentricity vector
% ecc_vect = 1/mu_E * cross(Y(:,4:6)', h_vect_t) - Y(:,1:3)'/vecnorm(Y(:,1:3)');  scritta così la cana
ecc_vect = 1/mu_E .* cross(v, h_vect) - r./vecnorm(r);
ecc_norm = vecnorm(ecc_vect);
e_h_prod = dot(ecc_vect, h_vect);

% radial and transversal velocity
ur = r./vecnorm(r);
uh = h_vect./h_norm;
ut = cross(uh, ur);

vr = dot(v, ur);
vt = dot(v, ut);


figure(2)
plot(T_ode, h_vect(1, :), 'b--', LineWidth=2);
hold on
plot(T_ode, h_vect(2, :), 'r--', LineWidth=2);
plot(T_ode, h_vect(3, :), 'g--', LineWidth=2);
plot(T_ode, h_norm, 'k--', LineWidth=2);
grid on
title('Angular momentum conservation');
xlabel('Time [s]');
ylabel('h_x, h_y, h_x, ||h|| [km^2/s]')

figure(3)
plot(T_ode, en);
title('Energy conservation');
xlabel('Time [s]');
ylabel('\epsilon [km^2/s^2]');
grid on

figure(4)
plot(T_ode, ecc_vect(1, :));
hold on
plot(T_ode, ecc_vect(2, :));
plot(T_ode, ecc_vect(3, :));
plot(T_ode, ecc_norm);
title('eccentricity')
legend('e_x', 'e_y', 'e_z', '||e||');
grid on

figure(5)
plot(T_ode, vecnorm(r))
title('position')
grid on
%axis equal

figure(6)
plot(T_ode, e_h_prod);
title('e*h')
grid on

figure
plot(T_ode, vr, 'r', LineWidth=2);
hold on
plot(T_ode, vt, 'b', LineWidth=2);
legend('v_r', 'v_{\theta}');
title('v_r and v_{\theta}');
grid on


%% EX1B_EARTH_2_BODY_PB
clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [6495; -970; -3622]; %[km]
v0 = [4.752; 2.130; 7.950]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
tspan = linspace(0, 4*T, 1000); % span for the ODE resolution

% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); % call this when you use the ode to set your own tollerances (adequate for the precision required for orbital porbelms

% perform the integration
[T_ode, Y] = ode113(@(t,y) ode_2pb(t, y, mu_E), tspan, y0, options);
                
                % choice of 113: more efficient than 45 with stingent tolleances
                % put the variables and call the system of equations
                % give the interval of time: t0 ti in which you want the evaluation tf
                % give the state vector 0
                % call the options about the tollerances
                % T = time, Y = matrix of the coordinatez x y z of the position at each time step asked

% Plot the results of the orbit
figure (1)
plot3(Y(:,1), Y(:,2), Y(:, 3), '-');
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit - Earth');
axis equal;
grid on;

% check the constants
r = Y(:, 1:3)';
v = Y(:, 4:6)';

% angular momentum
h_vect_0 = cross(r0, v0);
% h_vect_t = cross(Y(:,1:3)',Y(:,4:6)'); 
h_vect = cross(r, v);
h_norm = vecnorm(h_vect);

% energy
% en = dot(Y(:,4:6)', Y(:,4:6)')./2 - mu_E./vecnorm(Y(:,1:3)');
en = vecnorm(v).^2./2 - mu_E./vecnorm(r);

% eccentricity vector
% ecc_vect = 1/mu_E * cross(Y(:,4:6)', h_vect_t) - Y(:,1:3)'/vecnorm(Y(:,1:3)');  scritta così la cana
ecc_vect = 1/mu_E .* cross(v, h_vect) - r./vecnorm(r);
ecc_norm = vecnorm(ecc_vect);
e_h_prod = dot(ecc_vect, h_vect);

% radial and transversal velocity
ur = r./vecnorm(r);
uh = h_vect./h_norm;
ut = cross(uh, ur);

vr = dot(v, ur);
vt = dot(v, ut);


figure(2)
plot(T_ode, h_vect(1, :), 'b--', LineWidth=2);
hold on
plot(T_ode, h_vect(2, :), 'r--', LineWidth=2);
plot(T_ode, h_vect(3, :), 'g--', LineWidth=2);
plot(T_ode, h_norm, 'k--', LineWidth=2);
grid on
title('Angular momentum conservation');
xlabel('Time [s]');
ylabel('h_x, h_y, h_x, ||h|| [km^2/s]')

figure(3)
plot(T_ode, en);
title('Energy conservation');
xlabel('Time [s]');
ylabel('\epsilon [km^2/s^2]');
grid on

figure(4)
plot(T_ode, ecc_vect(1, :));
hold on
plot(T_ode, ecc_vect(2, :));
plot(T_ode, ecc_vect(3, :));
plot(T_ode, ecc_norm);
title('eccentricity')
legend('e_x', 'e_y', 'e_z', '||e||');
grid on

figure(5)
plot(T_ode, vecnorm(r))
title('position')
grid on
%axis equal

figure(6)
plot(T_ode, e_h_prod);
title('e*h')
grid on

figure
plot(T_ode, vr, 'r', LineWidth=2);
hold on
plot(T_ode, vt, 'b', LineWidth=2);
legend('v_r', 'v_{\theta}');
title('v_r and v_{\theta}');
grid on

%% EX2_PERTURBATED_2BP

%clear all
%close all
%clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]
R_E = astroConstants(23);
J_2 = astroConstants(9);

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [26578.137; 0; 0]; %[km]
v0 = [0; 2.221; 3.173]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
year = 365*24*60*60;
tspan_year = linspace(0, year, 1000); % span for the ODE resolution
nn = year/T;
NN = ceil(nn);
% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); % call this when you use the ode to set your own tollerances (adequate for the precision required for orbital porbelms

% perform the integration
[T_ode, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu_E, R_E, J_2), tspan_year, y0, options);

% T_mat = [];
% Y_mat = [];
% for i = 1:NN
%     t_span = linspace((i-1)*T, i*T, 1000);
%     [T_ode, Y_ode] = ode113(@(t,y) ode_pert_2pb(t, y, mu_E, R_E, J_2), t_span, y0, options);
%     T_mat = [T_mat; T_ode];
%     Y_mat = [Y_mat; Y_ode];
%     %i = i+1;
%     y0 = [Y_ode(end,:)'];
% end
%%
         
% Plot the results of the orbit
figure (1)
plot3(Y_mat(:,1), Y_mat(:,2), Y_mat(:, 3), '-');
xlabel('X [km]'); 
ylabel('Y [km]');
zlabel('Z [km]');
title('Perturbed Two-body problem orbit - Earth');
%axis equal;
grid on;

%%
% check the constants
r = Y(:, 1:3)';
v = Y(:, 4:6)';

% angular momentum
h_vect_0 = cross(r0, v0);
h_vect = cross(r, v);
h_norm = vecnorm(h_vect);

% energy
en = vecnorm(v).^2./2 - mu_E./vecnorm(r);

% eccentricity vector
ecc_vect = 1/mu_E .* cross(v, h_vect) - r./vecnorm(r);
ecc_norm = vecnorm(ecc_vect);
e_h_prod = dot(ecc_vect, h_vect);

% radial and transversal velocity
ur = r./vecnorm(r);
uh = h_vect./h_norm;
ut = cross(uh, ur);

vr = dot(v, ur);
vt = dot(v, ut);


figure(2)
plot(T_ode, h_vect(1, :), 'b--', LineWidth=2);
hold on
plot(T_ode, h_vect(2, :), 'r--', LineWidth=2);
plot(T_ode, h_vect(3, :), 'g--', LineWidth=2);
plot(T_ode, h_norm, 'k--', LineWidth=2);
grid on
title('Angular momentum conservation');
xlabel('Time [s]');
ylabel('h_x, h_y, h_x, ||h|| [km^2/s]')

figure(3)
plot(T_ode, en);
title('Energy conservation');
xlabel('Time [s]');
ylabel('\epsilon [km^2/s^2]');
grid on

figure(4)
plot(T_ode, ecc_vect(1, :));
hold on
plot(T_ode, ecc_vect(2, :));
plot(T_ode, ecc_vect(3, :));
plot(T_ode, ecc_norm);
title('eccentricity')
legend('e_x', 'e_y', 'e_z', '||e||');
grid on

figure(5)
plot(T_ode, vecnorm(r))
title('position')
grid on
%axis equal

figure(6)
plot(T_ode, e_h_prod);
title('e*h')
grid on

figure
plot(T_ode, vr, 'r', LineWidth=2);
hold on
plot(T_ode, vt, 'b', LineWidth=2);
legend('v_r', 'v_{\theta}');
title('v_r and v_{\theta}');
grid on

%% ES 1 - re-do

clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [26578.137; 0; 0]; %[km]
v0 = [0; 2.221; 3.173]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
%year = 365*24*60*60;
tspan = linspace(0, 2*T, 1000); % span for the ODE resolution

Orbit_Analysis(r0, v0, mu_E, tspan, 'non_perturbed', 'plot_orbit', 'more_plots');

%%
clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [6495; -970; -3622]; %[km]
v0 = [4.752; 2.130; 7.950]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
%year = 365*24*60*60;
tspan = linspace(0, 2*T, 1000); % span for the ODE resolution

Orbit_Analysis(r0, v0, mu_E, tspan, 'non_perturbed', 'plot_orbit', 'more_plots');

%% Perturbed 

clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]
R_E = astroConstants(23);
J_2 = astroConstants(9);

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [26578.137; 0; 0]; %[km]
v0 = [0; 2.221; 3.173]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
year = 365*24*60*60;
tspan = linspace(0, year, 10000); % span for the ODE resolution

Orbit_Analysis(r0, v0, mu_E, tspan, 'perturbed', 'plot_orbit','more_plots',J_2, R_E);

%%
%% Perturbed 

clear all
close all
clc

% Physiscal parameters
mu_E = astroConstants(13);     % [km^3/s^2]
R_E = astroConstants(23);
J_2 = astroConstants(9);

% Initial conditions 
% (position, velocity and state vector as column)
r0 = [6495; -970; -3622]; %[km]
v0 = [4.752; 2.130; 7.950]; %[km/s]
y0 = [r0; v0];

% Set time span
a = 1/(2/norm(r0) - dot(v0,v0)/mu_E); %[km]
T = 2*pi*sqrt(a^3/mu_E); %[s]
year = 365*24*60*60;
tspan = linspace(0, year, 10000); % span for the ODE resolution

Orbit_Analysis(r0, v0, mu_E, tspan, 'perturbed', 'plot_orbit','more_plots',J_2, R_E);

