function lambertOrbitPlotter(planet_dep, planet_arr, dep_arr_dates, window_dep, window_arr)

% Lambert transfer arc with the minimum deltaV given a launch window
% and an optional arrival window
%
% PROTOTYPE:
% lambertOrbitPlotter(planet_dep, planet_arr, dep_arr_dates)
%
% INPUT:
% planet_dep[1] Planet index of the departure planet [-]
% planet_arr[1] Planet index of the arrival planet  [-]
%               ...(1 through 8, starting from mercury) 
%
% dep_arr_dates[2x1] Two-element vector where the first is the departure
%                    date and the second is arrival date
%
%
% OUTPUT:
% A neat little plot of the lambert transfer arc
%
% CONTRIBUTORS:
% Mert Emrem
%
% VERSIONS
% 2024-11-20: 1.0
%

%% Orbits Plot

time_resolution = 1000;

% Earth and mars initial and final state vectors (r_vec, v_vec) in cartesian
% dep: departure, arr: arrival, opt: optimized value
% _d: departure planet, _a: arrival planet

mu_Sun = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]

ephemeris_d_dep = uplanet(dep_arr_dates(1), planet_dep); % Departure planet ephemeris at departure
ephemeris_a_dep = uplanet(dep_arr_dates(1), planet_arr); % Arrival planet ephemeris at departure
ephemeris_a_arr = uplanet(dep_arr_dates(2), planet_arr); % Arrival planet ephemeris at arrival
[r_d_dep, v_d_dep] = kep2car(ephemeris_d_dep, mu_Sun); % Cartesian state vector of departure planet at departure
[r_a_dep, v_a_dep] = kep2car(ephemeris_a_dep, mu_Sun); % Cartesian state vector of arrival planet at departure

[r_a_arr, ~] = kep2car(ephemeris_a_arr, mu_Sun);

ToF = (dep_arr_dates(2) - dep_arr_dates(1))*86400;

[~, ~, ~, ~, VI, ~, ~, ~] = ...
                lambertMR(r_d_dep, r_a_arr, ToF, mu_Sun, 0, 0, 0);

%Lambert Arc Plot
ToF_arc = (dep_arr_dates(2)-dep_arr_dates(1))*86400; % seconds
tspan = linspace( 0, ToF_arc, time_resolution );

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0_arc = [r_d_dep; VI']';

dy = @(t,y) ode_2body_pb(t, y, mu_Sun);
[ ~, Y_arc ] = ode113(dy, tspan, y0_arc, options);

%% Full orbits

T_dep = Eph2OrbitalPeriod(r_d_dep, v_d_dep);
T_arr = Eph2OrbitalPeriod(r_a_dep, v_a_dep);

% Time spans for full orbits
orbit_time_dep = linspace(0, T_dep, time_resolution);   % [s]
orbit_time_arr = linspace(0, T_arr, time_resolution);   % [s]

% Departure planet full orbit data
y0_d_full = [r_d_dep; v_d_dep]';
dy = @(t, y) ode_2body_pb(t, y, mu_Sun);
[~, Y_d_full] = ode113(dy, orbit_time_dep, y0_d_full, options);

% Departure planet motion over transfer time data
y0_d = [r_d_dep; v_d_dep]';
dy = @(t,y) ode_2body_pb(t, y, mu_Sun);
[ ~, Y_d ] = ode113(dy, tspan, y0_d, options);

% Arrival planet full orbit data
[r_a_orbit, v_a_orbit] = kep2car(uplanet(dep_arr_dates(1), planet_arr), mu_Sun);
y0_a_full = [r_a_orbit; v_a_orbit]';
[~, Y_a_full] = ode113(dy, orbit_time_arr, y0_a_full, options);

% Arrival planet motion over transfer time data
y0_a = [r_a_dep; v_a_dep]';
dy = @(t,y) ode_2body_pb(t, y, mu_Sun);
[ ~, Y_a ] = ode113(dy, tspan, y0_a, options);


if nargin == 5
    % Convert departure and arrival dates to MJD2000
    date_dep1 = date2mjd2000(window_dep(1,:));
    date_dep2 = date2mjd2000(window_dep(2,:));
    tspan_dep = linspace(0, date_dep2-date_dep1, 200);
    
    date_arr1 = date2mjd2000(window_arr(1,:));
    date_arr2 = date2mjd2000(window_arr(2,:));
    tspan_arr = linspace(0, date_arr2-date_arr1, 200);

    ephemeris_d_dep = uplanet(window_dep(1), planet_dep); % Departure planet ephemeris at departure
    ephemeris_a_arr = uplanet(window_arr(1), planet_arr); % Arrival planet ephemeris at arrival 
    [r_d_dep_window, v_d_dep_window] = kep2car(ephemeris_d_dep, mu_Sun); % Cartesian state vector of departure planet at first possible launch
    [r_a_arr_window, v_a_arr_window] = kep2car(ephemeris_a_arr, mu_Sun); % Cartesian state vector of arrival planet at first possible arrival

    % Departure planet motion over departure window data
    y0_d_dep = [r_d_dep_window; v_d_dep_window]';
    dy = @(t,y) ode_2body_pb(tspan_dep, y, mu_Sun);
    [ ~, Y_d_dep ] = ode113(dy, tspan_dep*86400, y0_d_dep, options);

    % Arrival planet motion over arrival window data
    y0_a_arr = [r_a_arr_window; v_a_arr_window]';
    dy = @(t,y) ode_2body_pb(tspan_arr, y, mu_Sun);
    [ ~, Y_a_arr ] = ode113(dy, tspan_arr*86400, y0_a_arr, options);

end

% 1 km = ... AU
km_to_au = 6.6845871226706E-9; 

f1 = figure;

f1.Position = [500 0 1200 500]; % Adjust figure window size

% Lambert Transfer Arc
plot3(Y_arc(:, 1) * km_to_au, Y_arc(:, 2) * km_to_au, Y_arc(:, 3) * km_to_au, ...
    'g', 'LineWidth', 5, 'Color', [0.1569    0.5804    0.2118]);
hold on;

% Departure planet's Full Orbit
plot3(Y_d_full(:, 1) * km_to_au, Y_d_full(:, 2) * km_to_au, Y_d_full(:, 3) * km_to_au, ...
    'r--', 'LineWidth', 3, 'Color', [0.8784    0.2196    0.3294]);

% Departure planet motion over transfer time plot
plot3(Y_d(:, 1) * km_to_au, Y_d(:, 2) * km_to_au, Y_d(:, 3) * km_to_au, ...
    'b', 'LineWidth', 3);

% Arrival planet's Full Orbit
plot3(Y_a_full(:, 1) * km_to_au, Y_a_full(:, 2) * km_to_au, Y_a_full(:, 3) * km_to_au, ...
    'b--', 'LineWidth', 3, 'Color', [0.3490    0.4667    0.8510]);

% Arrival planet motion over transfer time plot
plot3(Y_a(:, 1) * km_to_au, Y_a(:, 2) * km_to_au, Y_a(:, 3) * km_to_au, ...
    'r', 'LineWidth', 3);


if nargin == 5
    % Departure planet's motion over departure window plot
    p_e_departure = plot3(Y_d_dep(:, 1) * km_to_au, Y_d_dep(:, 2) * km_to_au, Y_d_dep(:, 3) * km_to_au, ...
        'b', 'LineWidth', 15);
    p_e_departure.Color(4) = 0.2;

    % Arrival planet motion over arrival window plot
    p_m_arrival = plot3(Y_a_arr(:, 1) * km_to_au, Y_a_arr(:, 2) * km_to_au, Y_a_arr(:, 3) * km_to_au, ...
        'r', 'LineWidth', 15);
    p_m_arrival.Color(4) = 0.2;
end

% Transfer Epoch Positions
plot3(r_d_dep(1) * km_to_au, r_d_dep(2) * km_to_au, r_d_dep(3) * km_to_au, ...
    'bo', 'MarkerFaceColor', 'b', 'MarkerSize',15);
plot3(r_a_arr(1) * km_to_au, r_a_arr(2) * km_to_au, r_a_arr(3) * km_to_au, ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize',12);

% Axis Labels and Styling
xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]');
legend('Transfer arc', 'Arrival planet during transfer', 'Arrival window', 'Departure planet during transfer', 'Departure window');
grid on;
axis equal;
title('Interplanetary Lambert Transfer');
hold off;


    function orbital_period = Eph2OrbitalPeriod(position, velocity)

        mu = astroConstants(4); % Sun's gravitational parameter [km^3/s^2]
    
        % Calculate semi-major axis (a)
        r = norm(position); % Distance from the Sun (m)
        v = norm(velocity); % Speed (m/s)
        specific_energy = 0.5 * v^2 - mu / r; % Orbital specific energy (J/kg)
        a = -mu / (2 * specific_energy); % Semi-major axis (m)

        % Compute the orbital period (T)
        orbital_period = 2 * pi * sqrt(a^3 / mu); % Orbital period (s)
    end

end

