function DeltaV_3dofs_Plotter(M, treshold, limit)


figure

[x, y, z] = meshgrid(1:size(M,1), 1:size(M,2), 1:size(M,3));  % Create the grid of x, y, and z

% Flatten the 3D arrays into 1D vectors for plotting
x = x(:);
y = y(:);
z = z(:);
values = M(:);   % Flatten BEH into a 1D vector for color mapping

valid_indices = values <= treshold;  % Logical array for valid points
x = x(valid_indices);        % Filtered x values
y = y(valid_indices);        % Filtered y values
z = z(valid_indices);        % Filtered z values
values = values(valid_indices);  % Filtered color values

% 3D Scatter plot
scatter3(x, y, z, 10, values, 'filled');  % 10 is the marker size, values control color
colorbar;                % Add a colorbar to interpret the color scale
clim([0 limit])
colormap parula;           
xlabel('t dep');
ylabel('t flyby');
zlabel('t arr');
axis equal;
title('DeltaV');
grid on;