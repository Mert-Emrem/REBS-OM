clear 
close all
clc

addpath 'MATLAB functions-20241213'\timeConversion\time\
addpath 'MATLAB functions-20241213'\


%% Interplanetary Explorer Mission

% constant
mu_sun = astroConstants(4);
mu_merc = astroConstants(11);
mu_mars = astroConstants(14);

R_mars = 3390; % mars radius [km]

%% Departure planet: Mercury

T_syn_1 = 100.8882;

dep_date_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
dep_date_max = date2mjd2000([2035, 1, 1, 0, 0, 0]);
% dep_date_max = dep_date_min + 10*T_syn_1; % ToF_max

dt = 20;

dep_window = dep_date_min: dt :dep_date_max;


[kep, ~] = uplanet_vec(dep_window, 1);

a_merc = kep(1,1);

[r_dep, v_dep] = kep2car_vec(a_merc, kep(2,:), kep(3,:),...
             kep(4,:), kep(5,:), kep(6,:), mu_sun);

%% Flyby planet: Mars

flyby_window = dep_window;

[kep, ~] = uplanet_vec(flyby_window, 4);

a_mars = kep(1,1);

[r_mars, v_mars] = kep2car_vec(a_mars, kep(2,:), kep(3,:),...
             kep(4,:), kep(5,:), kep(6,:), mu_sun);


%% Arrival asteroid N.40

T_syn_2 = 4.1902;

arr_time_min = date2mjd2000([2030, 1, 1, 0, 0, 0]);
arr_time_max = date2mjd2000([2035, 1, 1, 0, 0, 0]);

arr_dt = 100;

arr_window = arr_time_min: arr_dt: arr_time_max;

[kep, f, M] = ephAsteroids_vec(arr_window, 40);

a_harmonia = kep(1);
i = kep(3).*ones(1,length(f));
OM = kep(4).*ones(1, length(f));
om =  kep(5).*ones(1, length(f));

[r_harm, v_harm] = kep2car_vec(kep(1), kep(2), i,...
            OM, om, f, mu_sun);

%%

deltaV_Merc_Mars = ones(length(dep_window), length(flyby_window))*101;
Vinf_minus = ones(3, length(flyby_window));

for i = 1:length(dep_window)
    for j = 1:length(flyby_window)

        % Compute the time of flight
        ToF = (flyby_window(j) - dep_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, ~, VI, VF, ~, ~] = ...
                lambertMR(r_dep(:, i) , r_mars(:,j) , ...
                ToF, mu_sun, 0, 0, 0);

            Vinf_minus(:,j) = VF' - v_mars(:,j);

            % Compute Delta-V
            deltaV_1 = vecnorm(VI' - v_dep(:, i)); % Departure Delta-V
            deltaV_2 = vecnorm(VF' - v_mars(:,j)); % Arrival Delta-V
            deltaV_Merc_Mars(i, j) = deltaV_1+deltaV_2;

        end
    end
end

Merc_Mars_3d = repmat(deltaV_Merc_Mars, [1, 1, length(arr_window)]);

M = Merc_Mars_3d;

[x, y, z] = meshgrid(1:size(M,1), 1:size(M,2), 1:size(M,3));  % Create the grid of x, y, and z

% Flatten the 3D arrays into 1D vectors for plotting
x = x(:);
y = y(:);
z = z(:);
values = M(:);   % Flatten BEH into a 1D vector for color mapping

valid_indices = values <= 45;  % Logical array for valid points
x = x(valid_indices);        % Filtered x values
y = y(valid_indices);        % Filtered y values
z = z(valid_indices);        % Filtered z values
values = values(valid_indices);  % Filtered color values

% 3D Scatter plot
scatter3(x, y, z, 10, values, 'filled');  % 10 is the marker size, values control color
colorbar;                % Add a colorbar to interpret the color scale
clim([0 60])
colormap parula;            % Choose a colormap (e.g., 'jet', 'hot', 'parula')
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis equal;
title('3D Visualization of BEH with Color Representing 4th Dimension');
grid on;

porkchopPlotter(deltaV_Merc_Mars, flyby_window, dep_window)



deltaV_Mars_Harm = ones(length(flyby_window), length(arr_window))*101;
Vinf_plus = ones(3, length(flyby_window));

for i = 1:length(flyby_window)
    for j = 1:length(arr_window)

        % Compute the time of flight
        ToF = (arr_window(j) - flyby_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, ~, VI, VF, ~, ~] = ...
                lambertMR(r_mars(:, i) , r_harm(:,j) , ...
                ToF, mu_sun, 0, 0, 0);

            Vinf_plus(:,j) = VI' - v_mars(:,j);

            % Compute Delta-V
            deltaV_1 = vecnorm(VI' - v_mars(:,i)); % Departure Delta-V
            deltaV_2 = vecnorm(VF' - v_harm(:,j)); % Arrival Delta-V
            deltaV_Mars_Harm(i, j) = deltaV_1+deltaV_2;

        end
    end
end

%% 

Vinf_minus_val = vecnorm(Vinf_minus,2,1);
Vinf_plus_val = vecnorm(Vinf_plus,2,1);

delta_val = acos(dot(Vinf_minus, Vinf_plus)./(Vinf_minus_val.*Vinf_plus_val));

e = @(r_p, v_inf) 1+ (r_p.*(v_inf).^2)/mu_mars;

delta = @(r_p, v_inf) 2*asin(1./e(r_p, v_inf));

f_delta = @(r_p)...
            delta_val - (delta(r_p, Vinf_minus_val)/2 + delta(r_p, Vinf_plus_val)/2);

h_atm = 100000;

r_min = (R_mars + h_atm).*ones(1,length(flyby_window));
tol = 1e-14;
options = optimoptions('fsolve', 'TolFun', tol, 'TolX', tol);
r_p = fsolve(f_delta, r_min, options);

e_minus = e(r_p, Vinf_minus_val);
a_minus = r_p/(1-e_minus);

e_plus = e(r_p, Vinf_plus_val);
a_plus = r_p/(1-e_plus);

delta_v_flyby = Vinf_minus_val- Vinf_plus_val;

h_GA = r_p - R_mars;

v_p_minus = sqrt(Vinf_minus_val.^2 + 2*mu_mars./r_p);
v_p_plus = sqrt(Vinf_plus_val.^2 + 2*mu_mars./r_p);

delta_V_poweredFB = abs(v_p_plus - v_p_minus);

Mars_Harm_3d = repmat(deltaV_Mars_Harm, [1, 1, length(dep_window)]);

porkchopPlotter(deltaV_Mars_Harm, arr_window, dep_window)




%% plot
% 
% planet = 'Sun';
% opts.Units = 'km';
% opts.Position = [0, 0, 0];

% figure
% planet3D(planet, opts);
% view([54, 32])
% hold on
% grid on
% 
% x = r_dep(1,1);
% y = r_dep(2,1);
% z = r_dep(3,1);
% 
% % line=animatedline(x, y, z,'color', '#ffff00', 'LineWidth',2);
% line = animatedline;
% 
% for t = 1:size(r_dep,2)
% 
%         x = r_dep(1, t);
%         y = r_dep(2, t);
%         z = r_dep(3, t);
% 
% 
%         addpoints(line, x, y, z);
% 
%         drawnow;
% 
%         pause(0.01);
% 
% end

% plot3(r_dep(1,:), r_dep(2, :), r_dep(3, :))
% plot3(r_arr(1,:), r_arr(2, :), r_arr(3, :))
% plot3(r_f(1,:), r_f(2, :), r_f(3, :))
% legend('Sun','Mercury', 'Mars', 'Asteroid N40')




%%