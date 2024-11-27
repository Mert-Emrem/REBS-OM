function porkchopPlotter(deltaV_totals,tspan_arr,tspan_dep,fval_unc)
    
    t_offset = datenum('2000-01-01');
    [X, Y] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset); 
    f = figure;

    % Subplot 1: 2D Contour Plot
    subplot(1, 2, 1);
    deltaV_totals(deltaV_totals > 99) = NaN;
    contourf(X, Y, deltaV_totals', 100, 'LineStyle', 'none'); % Smooth 2D contours
    clim([fval_unc 100]);
    colorbar;
    xlabel('Departure Date');
    ylabel('Arrival Date');
    title('2D Porkchop Plot');
    hold on
    
    % Time-of-Flight lines

    resolution = 1000;
    ToF_matrix = zeros(resolution, resolution);
    
    for i = 1:resolution
        for j = 1:resolution
            ToF_matrix(i, j) = (tspan_arr(j) - tspan_dep(i)); % ToF in days
        end
    end
    
    % Time-of-Flight lines in days, showing every multiple of ToF_levels
    ToF_levels = 100;
    ToF_lines_min = round(min(ToF_matrix,[],"all")/ToF_levels)*ToF_levels;
    ToF_lines_max = round(max(ToF_matrix,[],"all")/ToF_levels)*ToF_levels;
    ToF_lines = ToF_lines_min:ToF_levels:ToF_lines_max; 
    contour(X, Y, ToF_matrix', ToF_lines, 'k', 'ShowText', 'on'); 

    % x-axis values only properly display this way
    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to x-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(45);

%%%%%%%%%%%%%%%%% 3D SURFACE PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1, 2, 2);
    surfc(X, Y, deltaV_totals', 'EdgeColor', 'none'); % 3D surface plot
    colormap('parula'); 
    zlim([fval_unc, 100]);
    clim([fval_unc 100]);
    xlabel('Departure Date');
    ylabel('Arrival Date');
    zlabel('Delta-V (km/s)');
    title('3D Porkchop Plot');
    
    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 6); % Set tick positions
    set(gca, 'XTick', xticks); % Apply ticks to x-axis
    xticklabels = datestr(xticks, 'yyyy-mmm-dd'); % Generate readable date labels
    set(gca, 'XTickLabel', xticklabels); % Set the date labels explicitly
    xtickangle(-45);
    
    datetick('y', 'yyyy-mmm-dd', 'keeplimits'); % Format y-axis as dates
    ytickangle(45);
    view(3); % Adjust to a 3D perspective
    
    f.Position = [500 500 1200 500]; % Adjust figure window size
end
