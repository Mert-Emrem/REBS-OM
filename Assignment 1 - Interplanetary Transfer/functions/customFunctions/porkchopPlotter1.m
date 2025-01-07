function porkchopPlotter1(deltaV_totals, tspan_arr, tspan_dep)

% Brief explanation:
% The function creates a porkchop plot for interplanetary trajectory analysis, 
% visualizing the delta-V required for transfers over a range of departure and 
% flyby dates. Both 2D and 3D visualizations are generated.

% 
% PROTOTYPE:
%  porkchopPlotter1(deltaV_totals, tspan_arr, tspan_dep)
% 
% INPUT:
%  deltaV_totals  [LxM]   Matrix of Delta-V values for the first leg, where
%                         L corresponds to departure dates, and M corresponds 
%                         to flyby dates [km/s].
%  tspan_arr      [1xM]   Vector of arrival dates (mjd2000).
%  tspan_dep      [1xL]   Vector of departure dates (mjd2000).
%
% OUTPUT:
%  Generates:
%  - A 2D contour plot showing Delta-V values and time-of-flight lines.
%  - A 3D surface plot illustrating Delta-V values as a function of departure 
%    and arrival dates.
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------



    % Set LaTeX as the default interpreter for text and axes
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    % Offset for MJD2000 time format conversion
    t_offset = datenum('2000-01-01');
    [X, Y] = meshgrid(tspan_dep + t_offset, tspan_arr + t_offset); 
    
    % Create 2D contour plot
    f1 = figure;
    min_deltaV = min(deltaV_totals(:));

    % Subplot 1: 2D Contour Plot
    contourf(X, Y, deltaV_totals', 100, 'LineStyle', 'none'); % Smooth 2D contours
    clim([min_deltaV 50]);
    colormap('parula');
    c = colorbar;
    c.Label.String = '$\Delta V$ [km/s]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 22;

    xlabel('Departure Date from Mercury', 'fontsize', 22, 'interpreter', 'latex');
    ylabel('Flyby Date at Mars', 'fontsize', 22, 'interpreter', 'latex');
    %title('$\Delta V$ of the First Leg of the Transfer', 'fontsize', 22, 'interpreter', 'latex');
    hold on;

    % Add Time-of-Flight (ToF) lines
    ToF_matrix = zeros(length(tspan_dep), length(tspan_arr));
    for i = 1:length(tspan_dep)
        for j = 1:length(tspan_arr)
            ToF_matrix(i, j) = (tspan_arr(j) - tspan_dep(i));
        end
    end
    ToF_levels = ceil(ToF_matrix(1, end) / 1000) * 75;
    ToF_lines_max = round(max(ToF_matrix, [], 'all') / ToF_levels) * ToF_levels;
    ToF_lines = 0:ToF_levels:ToF_lines_max; 
    [~, h] = contour(X, Y, ToF_matrix', ToF_lines, 'k', 'ShowText', 'on');
    set(findobj(h, 'Type', 'text'), 'Interpreter', 'latex', 'Color', 'k');

    % Format axis ticks with LaTeX-compatible date strings
    yticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 7); 
    set(gca, 'YTick', yticks, 'YTickLabel', latexifyDates(yticks)); 
    ytickangle(45);
    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 7); 
    set(gca, 'XTick', xticks, 'XTickLabel', latexifyDates(xticks)); 
    xtickangle(45);

    % Save 2D plot
    f1.Position = [100 100 800 800];

    hold off;

    % Create 3D surface plot
    f2 = figure;
    hold on;
    grid on
    deltaV_totals(deltaV_totals > 60) = NaN;
    surfc(X, Y, deltaV_totals', 'EdgeColor', 'none');
    colormap('parula');
    zlim([min_deltaV 50]);
    clim([min_deltaV 50]);

    xlabel('Departure Date from Mercury', 'fontsize', 22, 'interpreter', 'latex');
    ylabel('Flyby Date at Mars', 'fontsize', 22, 'interpreter', 'latex');
    zlabel('$\Delta V$ [km/s]', 'fontsize', 22, 'interpreter', 'latex');
    %title('$\Delta V$ of the First Leg of the Transfer', 'fontsize', 28, 'interpreter', 'latex');

    % Format Y-axis tick labels
    yticks = linspace(min(tspan_arr + t_offset), max(tspan_arr + t_offset), 7); 
    yticklabels = latexifyDates(yticks);
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'TickLabelInterpreter', 'latex');
    ytickangle(45);
    xticks = linspace(min(tspan_dep + t_offset), max(tspan_dep + t_offset), 7); 
    xticklabels = latexifyDates(xticks);
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex');
    xtickangle(45);

    view(3);
    f2.Position = [100 100 800 800];
    saveas(gcf, 'porkchop1_3d.png');
    saveas(gcf, 'porkchop1_3d', 'epsc');
end

% Helper function to convert dates to LaTeX-compatible strings
function latexDates = latexifyDates(dates)
    dateStrings = datestr(dates, 'yyyy-mmm-dd');
    dateStringsCell = cellstr(dateStrings);
    latexDates = cellfun(@(x) ['\texttt{' x '}'], dateStringsCell, 'UniformOutput', false);
end
