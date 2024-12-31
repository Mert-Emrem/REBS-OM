function [M1, M2, Vinf_minus, Vinf_plus] = LambertArcsDeltaV_calculator(dep_window, flyby_window, arr_window, FlagPorkchopPlot, FlagDeltaV_3dofs_Plot)


mu_sun = astroConstants(4);

%% Departure planet: Mercury

% Synodic period between Mercury-Mars, calculated using Estimation_of_syn_periods.m
%  T_syn_1 = 100.8882; [days]

% Obtain keplerian elements of mercury for each departure date [6 x L]
[dep_kep, ~] = uplanet_vec(dep_window, 1);

% Semimajor axis of mercury from keplerian elements
a_merc = dep_kep(1,1);

% Cartesian states of mercury for each departure date [(3,3) x L]
[r_dep, v_dep] = kep2car_vec(a_merc, dep_kep(2,1), dep_kep(3,1),...
             dep_kep(4,1), dep_kep(5,1), dep_kep(6,:), mu_sun);

%% Flyby planet: Mars

% Flyby window needs not be same as the departure window,
% but it allows for an intuitive visualization

% Same procedure as Mercury
[flyby_kep, ~] = uplanet_vec(flyby_window, 4);

a_mars = flyby_kep(1,1);

[r_mars, v_mars] = kep2car_vec(a_mars, flyby_kep(2,:), flyby_kep(3,:),...
             flyby_kep(4,:), flyby_kep(5,:), flyby_kep(6,:), mu_sun);

%% Arrival asteroid N.40 (Harmonia)

% Synodic period between Mars-Harmonia, calculated using Estimation_of_syn_periods.m
%  T_syn_2 = 4.1902;  [years]


% Obtain keplerian ephemerides of Harmonia
[arr_kep, f, ~] = ephAsteroids_vec(arr_window, 40);

a_harmonia = arr_kep(1);
i = arr_kep(3).*ones(1,length(f));
OM = arr_kep(4).*ones(1, length(f));
om =  arr_kep(5).*ones(1, length(f));

% Cartesian states of Harmonia for each arrival date [(3,3) x L]
[r_harm, v_harm] = kep2car_vec(arr_kep(1), arr_kep(2), i,...
            OM, om, f, mu_sun);

%% First leg of the transfer (Mercury - Mars)

% Preallocation of the matrices to be created
deltaV_Merc_Mars = ones(length(dep_window), length(flyby_window))*NaN;
Vinf_minus = ones(length(dep_window), length(flyby_window),3)*NaN;

for i = 1:length(dep_window)
    for j = 1:length(flyby_window)

        % Compute the time of flight
        ToF = (flyby_window(j) - dep_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (flyby before departure)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, err_1, VI, VF, ~, ~] = ...
                lambertMR(r_dep(:, i) , r_mars(:,j) , ...
                ToF, mu_sun, 0, 0, 0);


            % Compute Delta-V
            deltaV_1 = norm(VI' - v_dep(:, i)); % Departure Delta-V
            % deltaV_2 = norm(VF' - v_mars(:,j)); % Arrival Delta-V
            deltaV_2 = 0;

            if err_1 ~= 1 && err_1 ~= 3 && err_1 ~= 4
            deltaV_Merc_Mars(i, j) = deltaV_1+deltaV_2;
            Vinf_minus(i,j,:) = VF' - v_mars(:,j);
            end
        end
    end
end



%% Second leg of the transfer (Mars - Harmonia)

% Preallocation of the matrices to be created
deltaV_Mars_Harm = ones(length(flyby_window), length(arr_window))*NaN;
Vinf_plus = ones(length(flyby_window), length(arr_window),3)*NaN;

for i = 1:length(flyby_window)

    for j = 1:length(arr_window)

        % Compute the time of flight
        ToF = (arr_window(j) - flyby_window(i)) * 86400; % [s]
        
        % Ignore temporally impossible solutions (arrival at Harmonia before Mars flyby)
        if ToF > 0
            % Solve Lambert's problem
            [~, ~, ~, err_2, VI, VF, ~, ~] = ...
                lambertMR(r_mars(:, i) , r_harm(:,j) , ...
                ToF, mu_sun, 0, 0, 0);

            % Compute Delta-V
            deltaV_1 = norm(VI' - v_mars(:,i)); % Departure Delta-V
            deltaV_2 = norm(VF' - v_harm(:,j)); % Arrival Delta-V
            deltaV_1 = 0;

            if err_2 ~= 1 && err_2 ~= 3 && err_2 ~= 4
               deltaV_Mars_Harm(i, j) = deltaV_1+deltaV_2;
               Vinf_plus(i,j,:) = VI' - v_mars(:,i);
            end
        end
    end

end


%% Post-processing of the first and seconds legs, Gravity Assist

% Repeat DeltaV results of the first leg over every arrival date
Merc_Mars_3d = repmat(deltaV_Merc_Mars, [1, 1, length(arr_window)]);

% Repeat DeltaV results of the second leg over every departure date
Mars_Harm_3d = repmat(deltaV_Mars_Harm, [1, 1, length(dep_window)]);

% Reorient the DeltaV matrix of first leg such that the two are combinable
Merc_Mars_3d = permute(Merc_Mars_3d, [2, 3, 1]);

% Briefer names for ease of use
M1 = Mars_Harm_3d;
M2 = Merc_Mars_3d;

%%
if FlagPorkchopPlot
porkchopPlotter1(deltaV_Merc_Mars, flyby_window, dep_window)
porkchopPlotter2(deltaV_Mars_Harm, arr_window, dep_window)
end

%%
% Generate 2D and 3D porkchop plots of both legs of the transfer
if FlagDeltaV_3dofs_Plot
DeltaV_3dofs_Plotter(deltaV_Merc_Mars, deltaV_Mars_Harm, dep_window, flyby_window, arr_window)
end

end