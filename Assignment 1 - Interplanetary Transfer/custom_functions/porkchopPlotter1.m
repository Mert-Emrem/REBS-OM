function porkchopPlotter1(deltaV_totals, tspan_arr, tspan_dep)
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    t_offset = datenum('2000-01-01');
    [X, Y] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset); 
    f1 = figure;

    % Subplot 1: 2D Contour Plot
    contourf(X, Y, deltaV_totals', 100, 'LineStyle', 'none'); % Smooth 2D contours
    clim([0 120]);
    colormap('parula');
    c = colorbar;
    c.Label.String = '$\Delta V$ [km/s]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 14; % Set font size

    xlabel('Departure Date from Mercury', 'fontsize', 14, 'interpreter', 'latex');
    ylabel('Flyby Date at Mars', 'fontsize', 14, 'interpreter', 'latex');
    title('$\Delta V$ of the First Leg of the Transfer', 'fontsize', 14, 'interpreter', 'latex');
    hold on;

    % Time-of-Flight lines
    ToF_matrix = zeros(length(tspan_dep), length(tspan_arr));
    for i = 1:length(tspan_dep)
        for j = 1:length(tspan_arr)
            ToF_matrix(i, j) = (tspan_arr(j) - tspan_dep(i)); % ToF in days
        end
    end

    % Time-of-Flight lines in days, showing every multiple of ToF_levels
    ToF_levels = ceil(ToF_matrix(1,end)/1000)*75;
    ToF_lines_max = round(max(ToF_matrix, [], "all")/ToF_levels)*ToF_levels;
    ToF_lines = 0:ToF_levels:ToF_lines_max; 
    [~, h] = contour(X, Y, ToF_matrix', ToF_lines, 'k', 'ShowText', 'on');
    
    % Customize Time-of-Flight text labels
    set(findobj(h, 'Type', 'text'), 'Interpreter', 'latex', 'Color', 'k');

    % Format axis ticks with LaTeX-compatible date strings
    yticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 7); 
    set(gca, 'YTick', yticks, 'YTickLabel', latexifyDates(yticks)); 
    ytickangle(45);

    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 7); 
    set(gca, 'XTick', xticks, 'XTickLabel', latexifyDates(xticks)); 
    xtickangle(45);

    f1.Position = [100 100 800 800]; % Adjust figure window size
    saveas(gcf, 'porkchop1.png');
    saveas(gcf, 'porkchop1', 'epsc');
    hold off;

    % 3D Surface Plot
    f2 = figure;
    hold on;
    deltaV_totals(deltaV_totals > 50) = NaN;
    surfc(X, Y, deltaV_totals', 'EdgeColor', 'none'); % 3D surface plot
    colormap('parula'); 
    zlim([0 75]);
    clim([0 75]);

    xlabel('Departure Date from Mercury', 'fontsize', 14, 'interpreter', 'latex');
    ylabel('Flyby Date at Mars', 'fontsize', 14, 'interpreter', 'latex');
    zlabel('$\Delta V$ [km/s]', 'fontsize', 14, 'interpreter', 'latex');
    title('$\Delta V$ of the First Leg of the Transfer', 'fontsize', 14, 'interpreter', 'latex');

    % Format Y-axis tick labels with LaTeX-compatible date strings
    yticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 7); 
    yticklabels = latexifyDates(yticks); % Generate LaTeX-compatible labels
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'TickLabelInterpreter', 'latex');
    ytickangle(45);
    
    % Format X-axis tick labels with LaTeX-compatible date strings
    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 7); 
    xticklabels = latexifyDates(xticks); % Generate LaTeX-compatible labels
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex');
    xtickangle(45);

    view(3); % Adjust to a 3D perspective
    f2.Position = [100 100 800 800]; % Adjust figure window size
    saveas(gcf, 'porkchop1_3d.png');
    saveas(gcf, 'porkchop1_3d', 'epsc');
end

% Helper function to convert dates to LaTeX-compatible strings
function latexDates = latexifyDates(dates)
    % Convert numerical dates to formatted strings
    dateStrings = datestr(dates, 'yyyy-mmm-dd'); % Generate readable date strings
    dateStringsCell = cellstr(dateStrings); % Convert to cell array
    % Replace '-' with LaTeX-compatible dash for each string
    latexDates = cellfun(@(x) ['\texttt{' x '}'], dateStringsCell, 'UniformOutput', false);
end
