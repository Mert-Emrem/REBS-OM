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

[kep, f, ~] = ephAsteroids_vec(arr_window, 40);

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

                for k = 1:length(arr_window)

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
end

Merc_Mars_3d = repmat(deltaV_Merc_Mars, [1, 1, length(arr_window)]);
porkchopPlotter(deltaV_Merc_Mars, flyby_window, dep_window)