function Animated_Transfers_Plot(x)
% Function to plot animated interplanetary transfer arcs
% This function visualizes two transfer arcs, given the times of
% departure, flyby, and arrival.
% 
% PROTOTYPE:
%  Animated_Transfers_Plot(x)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: 
%                 x(1) = departure time (MJD2000)
%                 x(2) = flyby time (MJD2000)
%                 x(3) = arrival time (MJD2000)
% 
% OUTPUT:
%  Animated plot of three orbits and the transfer arcs
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Extract departure, flyby, and arrival times from the input array
t_d = x(1); % Departure time
t_f = x(2); % Flyby time
t_a = x(3); % Arrival time

% Compute the orbital elements and state vectors for the planets and the target asteroid
[kep_d, ksun] = uplanet(t_d, 1); % Orbital elements of Mercury at departure
[rr_d, v_merc] = kep2car(kep_d(1), kep_d(2), kep_d(3), kep_d(4), kep_d(5), kep_d(6), ksun);

[kep_f, ksun] = uplanet(t_f, 4); % Orbital elements of Mars at flyby
[rr_f, v_mars] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

[kep_a, ~, ~] = ephAsteroids(t_a, 40); % Orbital elements of asteroid at arrival
[rr_a, v_harm] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

[kep_f, ksun] = uplanet(t_d, 2); % Orbital elements of Venus at departure
[rr_d_ven, v_venus] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

[kep_f, ksun] = uplanet(t_d, 3); % Orbital elements of Earth at departure
[rr_d_earth, v_earth] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

% Solve Lambert's problem to determine transfer velocities
[~, ~, ~, ~, vt1_i, ~, ~, ~] = lambertMR(rr_d, rr_f, (t_f - t_d) * 24 * 3600, ksun, 0, 0, 0, 0);
[~, ~, ~, ~, vt2_i, ~, ~, ~] = lambertMR(rr_f, rr_a, (t_a - t_f) * 24 * 3600, ksun, 0, 0, 0, 0);

% Integrate the trajectories of the planets and asteroid over time
[~, r1] = TwoBodyPb([0 88*24*3600], [rr_d; v_merc], ksun); % Mercury
[~, r2] = TwoBodyPb([0 690*24*3600], [rr_f; v_mars], ksun); % Mars
[~, r3] = TwoBodyPb([0 1260*24*3600], [rr_a; v_harm], ksun); % Asteroid
[~, r4] = TwoBodyPb([0 366*24*3600], [rr_d_earth; v_earth], ksun); % Earth
[~, r5] = TwoBodyPb([0 226*24*3600], [rr_d_ven; v_venus], ksun); % Venus

% Convert trajectories to column vectors for plotting
r1 = r1';
r2 = r2';
r3 = r3';
r4 = r4';
r5 = r5';

[kep_f, ksun] = uplanet(t_d, 4); % Orbital elements of Mars at departure
[rr_d_mars, v_mars] = kep2car(kep_f(1), kep_f(2), kep_f(3), kep_f(4), kep_f(5), kep_f(6), ksun);

[kep_a, ~, ~] = ephAsteroids(t_d, 40); % Orbital elements of asteroid at departure
[rr_d_harm, v_harm] = kep2car(kep_a(1), kep_a(2), kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

% Define physical properties for objects
Mars.Radius = 3390;
Mars.name = 'Mars';

Mercury.Radius = 2439.7;
Mercury.name = 'Mercury';

Harmonia.Radius = 107.6;
Harmonia.name = 'Harmonia';

Sun.Radius = 700e+3;
Sun.name = 'Sun';

Venus.Radius = 6052;
Venus.name = 'Venus';

Earth.Radius = 6371;
Earth.name = 'Earth';

% Load the image
background = imread('MilkyWay.jpg');

%% Animated Plot:

figure

ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
imshow(background, 'Parent', ax1)

ax2 = axes('Position', [0 0 1 1]);
set(ax2, 'Color', 'none'); % Sfondo trasparente
hold(ax2, 'on');

view(ax2, 45, 25); % Set the viewing angle
axis(ax2, 'equal');
xlim(ax2,[-0.39e+9 +0.41e+9]);
ylim(ax2,[-0.39e+9 +0.41e+9]);
grid(ax2, 'on'); % Mostra la griglia
set(ax2, 'GridColor', 'w', 'GridAlpha', 0.25); % Colore della griglia bianco, semi-trasparente

% Add plot labels and title
xlabel(ax2,'X [Km]', 'fontsize', 14, 'interpreter', 'latex', 'Color', 'white');
ylabel(ax2,'Y [Km]', 'fontsize', 14, 'interpreter', 'latex', 'Color', 'white');
zlabel(ax2,'Z [Km]', 'fontsize', 14, 'interpreter', 'latex', 'Color', 'white');

% title
title(ax2,'Interplanetary Trajectory', 'fontsize', 20, 'interpreter', 'latex', 'Color', 'w');

% Plot the orbits of Mercury, Mars, and the asteroid
plot3(ax2,r1(1,:), r1(2, :), r1(3, :), 'LineStyle', '- -', 'Color', '#8fe866'); % Mercury
plot3(ax2,r2(1,:), r2(2, :), r2(3, :), 'LineStyle', '- -', 'Color', '#eb8552'); % Mars
plot3(ax2,r3(1,:), r3(2, :), r3(3, :), 'LineStyle', '- -', 'Color', '#bfb588'); % Asteroid
plot3(ax2,r4(1,:), r4(2, :), r4(3, :), 'LineStyle', '- -', 'Color', '#7295ed', 'LineWidth', 0.1); % Earth
plot3(ax2,r5(1,:), r5(2, :), r5(3, :), 'LineStyle', '- -', 'Color', '#f7bb72', 'LineWidth', 0.1); % Venus

scatter3(rr_d(1), rr_d(2), rr_d(3),40, 'r', 'filled')
scatter3(rr_f(1), rr_f(2), rr_f(3),40, 'b', 'filled')
scatter3(rr_a(1), rr_a(2), rr_a(3),40, 'm', 'filled')


% Add Lambert transfer arcs
line_y1 = animatedline(ax2,rr_d(1), rr_d(2), rr_d(3), 'color', '#0714fa', 'LineWidth', 2.5);
line_y2 = animatedline(ax2,rr_f(1), rr_f(2), rr_f(3), 'color', '#e607fa', 'LineWidth', 2.5);

% Add animated lines for the orbits
line1 = animatedline(ax2,rr_d(1), rr_d(2), rr_d(3), 'color', '#8fe866', 'LineWidth', 1.5);
line2 = animatedline(ax2,rr_d_mars(1), rr_d_mars(2), rr_d_mars(3), 'color', '#eb8552', 'LineWidth', 1.5);
line3 = animatedline(ax2,rr_d_harm(1), rr_d_harm(2), rr_d_harm(3), 'color', '#bfb588', 'LineWidth', 1.5);
line4 = animatedline(ax2,rr_d_earth(1), rr_d_earth(2), rr_d_earth(3), 'color', '#7295ed', 'LineWidth', 1);
line5 = animatedline(ax2,rr_d_ven(1), rr_d_ven(2), rr_d_ven(3), 'color', '#f7bb72', 'LineWidth', 1);

% Plot the Sun and initialize planet and Spacecraft animations
plotObjects(14, Sun, [0, 0, 0]);
[ObjVen, Ven_X, Ven_Y, Ven_Z] = Plot_Animated_Objects(1000, Venus, rr_d_ven);
[ObjMerc, Merc_X, Merc_Y, Merc_Z] = Plot_Animated_Objects(2000, Mercury, rr_d);
[ObjEarth, Ear_X, Ear_Y, Ear_Z] = Plot_Animated_Objects(1000, Earth, rr_d_earth);
[ObjMars, Mars_X, Mars_Y, Mars_Z] = Plot_Animated_Objects(1500, Mars, rr_d_mars);
[ObjHarm, Harm_X, Harm_Y, Harm_Z] = Plot_Animated_Objects(43000, Harmonia, rr_d_harm);
SpaceCraft = scatter3(ax2,rr_d(1), rr_d(2), rr_d(3), 30, 'r', '<', 'filled');

% Set real-time scaling factor
real_time_factor = 1e+10; % Adjust this to control animation speed

% Animate the transfer trajectory over time
ToF1 = (t_f - t_d) * 24 * 3600;
ToF2 = (t_a - t_f) * 24 * 3600;
dt = 86400/2; % 1/2 day in seconds
tspan = 0:dt:(ToF1 + ToF2);

% Initialize storage variables for trajectory points
r1 = [];
r2 = [];
r3 = [];
r4 = [];
r5 = [];
y = [];

% Add text for displaying the current time
depdate = mjd20002date(t_d);
textdep =sprintf('Time: $%02d:%02d:%02d$ $[h:m:s]$ \nDate: $%02d/%02d/%04d$ $[d/m/y]$', ...
                        depdate(4), depdate(5), ceil(depdate(6)), depdate(3), depdate(2), depdate(1));

% Update the date and time text color to light blue
time_text = text(-0.4e+9, 0.35e+9, textdep, 'FontSize', 16, 'Interpreter', 'latex', 'Color', 'w');

% Add legend to the plot
legend('Mercury Orbit', 'Mars Orbit', 'Harmonia Orbit', 'Earth Orbit', 'Venus Orbit',...
       'Departure', 'Flyby', 'Arrival','Lambert Arc (1)', 'Lambert Arc (2)', ...
       'fontsize', 15, 'interpreter', 'latex', 'TextColor', 'k');

% Pause for a moment before starting the animation
pause(3.5);

% Loop through time steps to animate the transfers and object orbits
for jj = 1:length(tspan) - 1
    % Determine the transfer arc (first or second)
    if tspan(jj) <= ToF1
        [~, y_temp] = TwoBodyPb([0 tspan(jj+1)], [rr_d; vt1_i'], ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y1, y(1,1), y(1,2), y(1,3));
    else
        [~, y_temp] = TwoBodyPb([ToF1 tspan(jj+1)], [rr_f; vt2_i'], ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y2, y(1,1), y(1,2), y(1,3));
    end


    % Update object positions for animation
    [~, rtemp1] = TwoBodyPb([0 tspan(jj+1)], [rr_d; v_merc], ksun);
    r1 = [rtemp1(end, 1:3); r1];
    [~, rtemp2] = TwoBodyPb([0 tspan(jj+1)], [rr_d_mars; v_mars], ksun);
    r2 = [rtemp2(end, 1:3); r2];
    [~, rtemp3] = TwoBodyPb([0 tspan(jj+1)], [rr_d_harm; v_harm], ksun);
    r3 = [rtemp3(end, 1:3); r3];
    [~, rtemp4] = TwoBodyPb([0 tspan(jj+1)], [rr_d_earth; v_earth], ksun);
    r4 = [rtemp4(end, 1:3); r4];
    [~, rtemp5] = TwoBodyPb([0 tspan(jj+1)], [rr_d_ven; v_venus], ksun);
    r5 = [rtemp5(end, 1:3); r5];

    % Update positions for animated objects
    set(SpaceCraft, 'XData', y(1,1), 'YData', y(1,2), 'ZData', y(1, 3));
    set(ObjMerc, 'XData', r1(1,1) + Merc_X, 'YData', r1(1,2) + Merc_Y, 'ZData', r1(1,3) + Merc_Z);
    set(ObjMars, 'XData', r2(1,1) + Mars_X, 'YData', r2(1,2) + Mars_Y, 'ZData', r2(1,3) + Mars_Z);
    set(ObjHarm, 'XData', r3(1,1) + Harm_X, 'YData', r3(1,2) + Harm_Y, 'ZData', r3(1,3) + Harm_Z);
    set(ObjEarth, 'XData', r4(1,1) + Ear_X, 'YData', r4(1,2) + Ear_Y, 'ZData', r4(1,3) + Ear_Z);
    set(ObjVen, 'XData', r5(1,1) + Ven_X, 'YData', r5(1,2) + Ven_Y, 'ZData', r5(1,3) + Ven_Z);
    
    addpoints(line1, r1(1,1), r1(1,2), r1(1,3));
    addpoints(line2, r2(1,1), r2(1,2), r2(1,3));
    addpoints(line3, r3(1,1), r3(1,2), r3(1,3));
    addpoints(line4, r4(1,1), r4(1,2), r4(1,3));
    addpoints(line5, r5(1,1), r5(1,2), r5(1,3));

    % Update the displayed time
    current_time_mjd2000 = t_d + tspan(jj) / (24 * 3600);
    time = mjd20002date(current_time_mjd2000);

    % Format the time string with leading zeros
    time_str =sprintf('Time: $%02d:%02d:%02d$ $[h:m:s]$ \nDate:  $%02d/%02d/%04d$ $[d/m/y]$', ...
                        time(4), time(5), ceil(time(6)), time(3), time(2), time(1));

    % Update time text
    set(time_text, 'String', time_str);

    % Adjust real-time scaling
    pause(dt / real_time_factor);

    drawnow;
end


hold off;
end
