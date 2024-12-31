function Total_DeltaV_intersections = DeltaV_3dofs_Plotter(deltaV_totals_1, deltaV_totals_2, tspan_dep, tspan_flyby, tspan_arr, Vinf_minus, Vinf_plus, data)
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

    delta_limit_leg_1 = 30;
    delta_limit_leg_2 = 10;

    t_offset = datenum('2000-01-01');
    [X, Y] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset); 

    % Limit values for visualization
    deltaV_totals_viable_1 = deltaV_totals_1;
    deltaV_totals_viable_1(deltaV_totals_1 > delta_limit_leg_1) = NaN; 

    % Limit values for visualization
    deltaV_totals_viable_2 = deltaV_totals_2;
    deltaV_totals_viable_2(deltaV_totals_2 > delta_limit_leg_2) = NaN; 

    % Create 3D Visualization for the first leg
    depth = tspan_arr + t_offset; 
    Z3D = repmat(deltaV_totals_viable_1, [1, 1, numel(depth)]); 

    % Create 3D matrix for the second leg
    dep_window = length(tspan_dep); % Number of layers along the X-axis
    Z3D_2 = permute(repmat(deltaV_totals_viable_2, [1, 1, dep_window]), [3, 1, 2]); % Repeat and permute

    %% First leg
    % Visualization
    f0 = figure;
    hold on;
    f0.Position = [300 100 800 600]; % Adjust figure window size
    hold on;

    % Add 2D porkchop plot for deltaV_totals_1 (constant Z at base)
    baseZ = depth(1); % Fix Z at the bottom
    surf(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9); % Project deltaV_totals_1
    c = colorbar; 
    c.Label.String = '$\Delta V \, \mathrm{[km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    contour3(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        50, 'k');
    
    % Visualization of deltaV_totals_1 (X-Y slices)
    [X3D_1, Y3D_1, Z3DDepth_1] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset, tspan_arr + t_offset);
    hSlices1 = slice(X3D_1, Y3D_1, Z3DDepth_1, Z3D, [], [], tspan_arr + t_offset); % Slices along Z
    shading interp; % Smooth shading
    alpha(hSlices1, 0.6); % Transparency
    clim([0 100]);

       % Set axis limits
    xlim([min(tspan_flyby + t_offset), max(tspan_flyby + t_offset)]); % X-axis limits
    ylim([min(tspan_dep + t_offset), max(tspan_dep + t_offset)]); % Y-axis limits
    zlim([min(tspan_arr + t_offset), max(tspan_arr + t_offset)]); % Z-axis limits
    
    % Set X-axis ticks and labels (Flyby dates)
    xticks = linspace(min(tspan_flyby + t_offset), max(tspan_flyby + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to X-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(45); % Rotate labels for better readability
    
    % Set Y-axis ticks and labels (Departure dates)
    yticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'YTick', yticks); % Apply ticks to Y-axis
    yticklabels = datestr(yticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'YTickLabel', yticklabels); % Set the date labels explicitly
    ytickangle(20); % Rotate labels for better readability
    
    % Set Z-axis ticks and labels (Arrival dates)
    zticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 6); % Set tick positions
    set(gca, 'ZTick', zticks); % Apply ticks to Z-axis
    zticklabels = datestr(zticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'ZTickLabel', zticklabels); % Set the date labels explicitly

    grid on;

    view(3); % Set 3D perspective
    saveas(gcf, 'deltaV_totals_1.png');
    saveas(gcf, 'deltaV_totals_1', 'epsc');

    %% Second leg
        
    f1 = figure;
    hold on;
    f1.Position = [300 100 800 600]; % Adjust window size
    % Generate 2D grids for X (departure dates) and Z (arrival dates)
    [X_grid, Z_grid] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset);
    
    % Add 2D porkchop plot for deltaV_totals_2 (constant Y at base)
    baseY = tspan_flyby(end) + t_offset; % Fix Y at its first value
    surf(X_grid, baseY * ones(size(deltaV_totals_2')), Z_grid, deltaV_totals_2', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9); % Project the viable data at constant Y
    contour3(X_grid, baseY * ones(size(deltaV_totals_2)), Z_grid, deltaV_totals_2, ...
        50, 'k'); % Add contour lines for clarity

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
    c = colorbar; % Move to the right
    c.Label.String = '$\Delta V \, \mathrm{[km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';

    clim([0 50]);
    view(3); % 3D view
    hold off;
    % Set axis limits
    xlim([min(tspan_flyby + t_offset), max(tspan_flyby + t_offset)]); % X-axis limits
    ylim([min(tspan_dep + t_offset), max(tspan_dep + t_offset)]); % Y-axis limits
    zlim([min(tspan_arr + t_offset), max(tspan_arr + t_offset)]); % Z-axis limits
    
    % Set X-axis ticks and labels (Flyby dates)
    xticks = linspace(min(tspan_flyby + t_offset), max(tspan_flyby + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to X-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(45); % Rotate labels for better readability
    
    % Set Y-axis ticks and labels (Departure dates)
    yticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'YTick', yticks); % Apply ticks to Y-axis
    yticklabels = datestr(yticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'YTickLabel', yticklabels); % Set the date labels explicitly
    ytickangle(20); % Rotate labels for better readability
    
    % Set Z-axis ticks and labels (Arrival dates)
    zticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 6); % Set tick positions
    set(gca, 'ZTick', zticks); % Apply ticks to Z-axis
    zticklabels = datestr(zticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'ZTickLabel', zticklabels); % Set the date labels explicitly
    grid on;
    view(3); % Set 3D perspective
    saveas(gcf, 'deltaV_totals_2.png');
    saveas(gcf, 'deltaV_totals_2', 'epsc');
    
    %% Combined

    f2 = figure;
    hold on;
    f2.Position = [300 100 800 600]; % Adjust window size
    hold on;
    % Add 2D porkchop plot for deltaV_totals_1 (constant Z at base)
    baseZ = depth(1); % Fix Z at the bottom
    surf(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        'EdgeColor', 'none', 'FaceAlpha', 1); % Project deltaV_totals_1
    colormap(gca, parula);
    c = colorbar; % Move to the right
    c.Label.String = '$\Delta V \, \mathrm{[km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';

    contour3(X, Y, baseZ * ones(size(deltaV_totals_1)), deltaV_totals_1, ...
        50, 'k');

    % Visualization of deltaV_totals_1 (X-Y slices)
    [X3D_1, Y3D_1, Z3DDepth_1] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset, tspan_arr + t_offset);
    hSlices1 = slice(X3D_1, Y3D_1, Z3DDepth_1, Z3D, [], [], tspan_arr + t_offset); % Slices along Z
    shading interp; % Smooth shading
    alpha(hSlices1, 0.6); % Transparency

    % Generate 2D grids for X (departure dates) and Z (arrival dates)
    [X_grid, Z_grid] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset);
    
    % Add 2D porkchop plot for deltaV_totals_2 (constant Y at base)
    baseY = tspan_flyby(end) + t_offset; % Fix Y at its first value
    surf(X_grid, baseY * ones(size(deltaV_totals_2')), Z_grid, deltaV_totals_2', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9); % Project the viable data at constant Y
    contour3(X_grid, baseY * ones(size(deltaV_totals_2)), Z_grid, deltaV_totals_2, ...
        50, 'k'); % Add contour lines for clarity

    % Visualization of deltaV_totals_2 (Y-Z slices)
    [X3D_2, Y3D_2, Z3DDepth_2] = meshgrid(tspan_dep + t_offset, tspan_flyby + t_offset, tspan_arr + t_offset);
    hSlices2 = slice(X3D_2, Y3D_2, Z3DDepth_2, Z3D_2, tspan_dep + t_offset, [], []); % Slices along X
    shading interp; % Smooth shading
    alpha(hSlices2, 0.25); % Transparency

    % Set axis limits
    xlim([min(tspan_flyby + t_offset), max(tspan_flyby + t_offset)]); % X-axis limits
    ylim([min(tspan_dep + t_offset), max(tspan_dep + t_offset)]); % Y-axis limits
    zlim([min(tspan_arr + t_offset), max(tspan_arr + t_offset)]); % Z-axis limits
    
    % Set X-axis ticks and labels (Flyby dates)
    xticks = linspace(min(tspan_flyby + t_offset), max(tspan_flyby + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to X-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(45); % Rotate labels for better readability
    
    % Set Y-axis ticks and labels (Departure dates)
    yticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'YTick', yticks); % Apply ticks to Y-axis
    yticklabels = datestr(yticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'YTickLabel', yticklabels); % Set the date labels explicitly
    ytickangle(20); % Rotate labels for better readability
    
    % Set Z-axis ticks and labels (Arrival dates)
    zticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 6); % Set tick positions
    set(gca, 'ZTick', zticks); % Apply ticks to Z-axis
    zticklabels = datestr(zticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'ZTickLabel', zticklabels); % Set the date labels explicitly


    % Labels and settings
    ylabel('Departure date from Mercury');
    xlabel('Flyby date at Mars');
    zlabel('Arrival date at A40');
    colormap(gca, parula); 

    clim([0 50])
    view(3); % 3D view
    hold off;

    saveas(gcf, 'deltaV_totals_combined.png');
    saveas(gcf, 'deltaV_totals_combined', 'epsc');


    %% Intersection
    % Find the intersection points directly from Z3D and Z3D_2
    intersection_mask = ~isnan(Z3D) & ~isnan(Z3D_2); % Logical mask for valid intersections
    
    % Extract coordinates and values at intersection points
    X_intersect = X3D_1(intersection_mask);
    Y_intersect = Y3D_1(intersection_mask);
    Z_intersect = Z3DDepth_1(intersection_mask);
    
    % Extract deltaV values at intersection points
    deltaV_intersect_1 = Z3D(intersection_mask);
    deltaV_intersect_2 = Z3D_2(intersection_mask);

    % Initialize arrays
    Delta_GA_intersections = NaN(size(X_intersect));
    Total_DeltaV_intersections = NaN(size(X_intersect));
    

    %%
    % Loop over intersection points
    for idx = 1:numel(X_intersect)
        fprintf('Current index: %d\n',idx);
        % Find the associated indices in the grid
        [dep_idx, flyby_idx, arr_idx] = ind2sub(size(Z3D), find(intersection_mask));
    
        % Extract velocity vectors
        vinf_m = squeeze(Vinf_minus(dep_idx(idx), flyby_idx(idx), :)); % First leg outgoing
        vinf_p = squeeze(Vinf_plus(flyby_idx(idx), arr_idx(idx), :)); % Second leg incoming
    
        % Perform gravity assist calculation
        [dvp, ~, rp] = PowerGravityAssist(vinf_m, vinf_p, data.Mars.Radius, data.Mars.h_atm, data.Mars.mu, 0);
    
        % Store gravity assist Delta-V if valid
        if rp > (data.Mars.Radius + data.Mars.h_atm)
            Delta_GA_intersections(idx) = dvp;
            Total_DeltaV_intersections(idx) = dvp + deltaV_intersect_1(idx) + deltaV_intersect_2(idx);
        end
    end
        

    % Plot intersection points in 3D
    f3 = figure;
    hold on;
    f3.Position = [300 100 800 600]; % Adjust window size
    scatter3(X_intersect, Y_intersect, Z_intersect, 50, Total_DeltaV_intersections, 'filled'); 
    % Find the index of the minimum Delta-V
    [min_DeltaV, min_idx] = min(Total_DeltaV_intersections);
    hold on;
    % Highlight the minimum point
    scatter3(X_intersect(min_idx), Y_intersect(min_idx), Z_intersect(min_idx), 100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    %%
    % Add lines connecting the minimum point to the axes
    line([X_intersect(min_idx), X_intersect(min_idx)], [Y_intersect(min_idx), Y_intersect(min_idx)], [0, Z_intersect(min_idx)], ...
        'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5); % Line to Z-axis
    line([X_intersect(min_idx), X_intersect(min_idx)], [0, Y_intersect(min_idx)], [Z_intersect(min_idx), Z_intersect(min_idx)], ...
        'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5); % Line to Y-axis
    line([0, X_intersect(min_idx)], [Y_intersect(min_idx), Y_intersect(min_idx)], [Z_intersect(min_idx), Z_intersect(min_idx)], ...
        'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5); % Line to X-axis
    
        % Calculate the range of X and Y coordinates
        x_range = max(tspan_flyby) - min(tspan_flyby);
        y_range = max(tspan_dep) - min(tspan_dep);
        z_range = max(tspan_arr) - min(tspan_arr);
        
    % Define the percentage offsets (between 0 and 1)
    percent_offset_x = 0.05; % 5% of the X range (adjust as needed)
    percent_offset_y = 0.05; % 5% of the Y range (adjust as needed)
    percent_offset_z = 0.05; % 5% of the Y range (adjust as needed)
    
    % Calculate the offset values based on the range and percentage
    offset_x = percent_offset_x * x_range;
    offset_y = percent_offset_y * y_range;
    offset_z = percent_offset_z * z_range;

    text(X_intersect(min_idx)+offset_x, Y_intersect(min_idx)+offset_y, Z_intersect(min_idx)+offset_z, ...
    sprintf('$$\\Delta V$$ = %.2f km/s', min_DeltaV), ...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    %%
    colormap(gca, parula);
    c = colorbar; % Move to the right
    c.Label.String = '$\Delta V \, \mathrm{[km/s]}$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    
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
    xtickangle(45); % Rotate labels for better readability
    
    % Set Y-axis ticks and labels (Departure dates)
    yticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'YTick', yticks); % Apply ticks to Y-axis
    yticklabels = datestr(yticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'YTickLabel', yticklabels); % Set the date labels explicitly
    ytickangle(20); % Rotate labels for better readability
    
    % Set Z-axis ticks and labels (Arrival dates)
    zticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 6); % Set tick positions
    set(gca, 'ZTick', zticks); % Apply ticks to Z-axis
    zticklabels = datestr(zticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'ZTickLabel', zticklabels); % Set the date labels explicitly

    grid on;

end
