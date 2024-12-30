function total_delta_V_intersect = DeltaV_3dofs_Plotter(deltaV_totals_1, deltaV_totals_2, tspan_dep, tspan_flyby, tspan_arr)
    % Function to compute the Powered Gravity Assist.
    % 
    % INPUT:
    %  deltaV_totals_1 [LxM] Matrix of deltaV amounts where L is the size of
    %                      departure dates and M is flyby dates [km/s]
    %  deltaV_totals_2 [MxN] Matrix of deltaV amounts where M is the size of
    %                      departure flyby and N is arrival dates [km/s]
    %  tspan_dep   [1xL]     Departure date vector mjd2000
    %  tspan_flyby [1xM]     Flyby date vector in mjd2000
    %  tspan_arr   [1xN]     Arrival date vector in mjd2000
    
    % OUTPUT:
    % 3D porkchop plot visualization

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    delta_limit_leg_1 = 20;
    delta_limit_leg_2 = 5;

    t_offset = datenum('2000-01-01');
    [X, Y] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset); 

    % Limit values for visualization
    deltaV_totals_viable_1 = deltaV_totals_1;
    deltaV_totals_viable_1(deltaV_totals_1 > delta_limit_leg_1) = NaN; 

    % Limit values for visualization
    deltaV_totals_viable_2 = deltaV_totals_2;
    deltaV_totals_viable_2(deltaV_totals_2 > delta_limit_leg_2) = NaN; 

    % Time-of-Flight lines
    ToF_matrix = zeros(length(tspan_dep), length(tspan_flyby));
    for i = 1:length(tspan_dep)
        for j = 1:length(tspan_flyby)
            ToF_matrix(i, j) = (tspan_flyby(j) - tspan_dep(i)); % ToF in days
        end
    end
        
    % Time-of-Flight levels
    ToF_levels = ceil(ToF_matrix(1,end)/1000)*75;
    ToF_lines_max = round(max(ToF_matrix,[],"all")/ToF_levels)*ToF_levels;
    ToF_lines = 0:ToF_levels:ToF_lines_max; 

    % Create 3D Visualization for the first leg
    depth = tspan_arr + t_offset; 
    Z3D = repmat(deltaV_totals_viable_1, [1, 1, numel(depth)]); 

    % Create 3D matrix for the second leg
    dep_window = length(tspan_dep); % Number of layers along the X-axis
    Z3D_2 = permute(repmat(deltaV_totals_viable_2, [1, 1, dep_window]), [3, 1, 2]); % Repeat and permute

    % Visualization
    figure;
    hold on;

    % Add 2D porkchop plot for deltaV_totals_1 (constant Z at base)
    baseZ = depth(1); % Fix Z at the bottom
    surf(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9); % Project deltaV_totals_1
    colormap(gca, gray);
    colorbar;
    contour3(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        50, 'k');
    clim([0 50]);
    
    % Generate 2D grids for X (departure dates) and Z (arrival dates)
    [X_grid, Z_grid] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset);
    
    % Add 2D porkchop plot for deltaV_totals_2 (constant Y at base)
    baseY = tspan_flyby(end) + t_offset; % Fix Y at its first value
    surf(X_grid, baseY * ones(size(deltaV_totals_2')), Z_grid, deltaV_totals_2', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9); % Project the viable data at constant Y
    colorbar; % Add colorbar for the base
    contour3(X_grid, baseY * ones(size(deltaV_totals_2)), Z_grid, deltaV_totals_2, ...
        50, 'k'); % Add contour lines for clarity
    clim([0 150]);

    % Visualization of deltaV_totals_1 (X-Y slices)
    [X3D_1, Y3D_1, Z3DDepth_1] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset, tspan_arr + t_offset);
    hSlices1 = slice(X3D_1, Y3D_1, Z3DDepth_1, Z3D, [], [], tspan_arr + t_offset); % Slices along Z
    shading interp; % Smooth shading
    alpha(hSlices1, 0.6); % Transparency
    
    % Visualization of deltaV_totals_2 (Y-Z slices)
    [X3D_2, Y3D_2, Z3DDepth_2] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset, tspan_arr + t_offset);
    hSlices2 = slice(X3D_2, Y3D_2, Z3DDepth_2, Z3D_2, tspan_dep + t_offset, [], []); % Slices along X
    shading interp; % Smooth shading
    alpha(hSlices2, 0.25); % Transparency
    
    % Labels and settings
    ylabel('Departure date from Mercury');
    xlabel('Flyby date at Mars');
    zlabel('Arrival date at A40');
    colormap(gca, parula); 
    colorbar; 
    clim([0 50]);
    c.Label.String = '$\Delta V$ [km/s]';
    c.Label.Interpreter = 'latex';
    view(3); % 3D view
    hold off;

    % Find the intersection points directly from Z3D and Z3D_2
    intersection_mask = ~isnan(Z3D) & ~isnan(Z3D_2); % Logical mask for valid intersections
    
    % Extract coordinates and values at intersection points
    X_intersect = X3D_1(intersection_mask);
    Y_intersect = Y3D_1(intersection_mask);
    Z_intersect = Z3DDepth_1(intersection_mask);
    
    % Extract deltaV values at intersection points
    deltaV_intersect_1 = Z3D(intersection_mask);
    deltaV_intersect_2 = Z3D_2(intersection_mask);
    total_delta_V_intersect = deltaV_intersect_1 + deltaV_intersect_2;
    
    % Plot intersection points in 3D
    figure;
    scatter3(X_intersect, Y_intersect, Z_intersect, 50, total_delta_V_intersect, 'filled'); % Color by deltaV from leg 1
    colorbar; % Add colorbar for reference
    colormap(parula); % Color map for deltaV visualization
    xlabel('Flyby date at Mars');
    ylabel('Departure date from Mercury');
    zlabel('Arrival date at A40');
    title('Intersection of Viable DeltaV Points');
    view(3); % 3D view
    % Set axis limits
    xlim([min(tspan_flyby + t_offset), max(tspan_flyby + t_offset)]); % X-axis limits
    ylim([min(tspan_dep + t_offset), max(tspan_dep + t_offset)]); % Y-axis limits
    zlim([min(tspan_arr + t_offset), max(tspan_arr + t_offset)]); % Z-axis limits
    
    % Set X-axis ticks and labels (Flyby dates)
    xticks = linspace(min(tspan_flyby + t_offset), max(tspan_flyby + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to X-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(-45); % Rotate labels for better readability
    
    % Set Y-axis ticks and labels (Departure dates)
    yticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'YTick', yticks); % Apply ticks to Y-axis
    yticklabels = datestr(yticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'YTickLabel', yticklabels); % Set the date labels explicitly
    ytickangle(45); % Rotate labels for better readability
    
    % Set Z-axis ticks and labels (Arrival dates)
    zticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 6); % Set tick positions
    set(gca, 'ZTick', zticks); % Apply ticks to Z-axis
    zticklabels = datestr(zticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'ZTickLabel', zticklabels); % Set the date labels explicitly

    grid on;



end
