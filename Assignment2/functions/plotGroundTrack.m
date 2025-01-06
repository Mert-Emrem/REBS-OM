function [] = plotGroundTrack(axes, lon, lat, t)
% -------------------------------------------------------------------------
% [] = plotGroundTrack(lon, lat, t)
%
% Function plotGroundTrack plots the groundtrack of a given orbit. It
% should be used after function groundTrack, wich give as output the input
% of this function. 
% Notice: upload the image of the planet in the same folder of this
% function.
%
% INPUT:
% - axis      [gca]     [-]          gca - to make it feasible also for subplots
% - lon       [nx1]     [deg]        longitude
% - lat       [nx1]     [deg]        latitude
% - t         [nx1]     [s]          column vector of time in witch to compute the solution
%        
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% Last update: 30/12/2024
% -------------------------------------------------------------------------

% Upload the image of Mercury
img = imread('MercuryTexture.jpg');

%[rows, cols, ~] = size(img);

% Create the plot
%figure;
ax = axes;

%Show the image as background and make it cover the span of the axis
imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', ax);
set(ax, 'YDir', 'normal');
hold on;
axis on;
xlabel('Longitude [deg]', 'Interpreter','latex', FontSize=22);
ylabel('Latitude [deg]', 'Interpreter','latex', FontSize=22);


% limits for the axis
xlim([-180 180]);
ylim([-90 90]);

% Greenwich and Equator
plot([0 0], ylim, 'k', 'LineWidth', 1); % Greenwich Meridian of the planet
plot(xlim, [0 0], 'k--', 'LineWidth', 1); % Equator

scatter(lon, lat, 5, t, 'filled');
colorbar
%plot(lon, lat, 'Color', 'y', LineWidth=1);
plot(lon(1), lat(1), 'or', LineWidth=2);
plot(lon(end), lat(end), 'ob', LineWidth=2);
ylabel(colorbar,'Time [s]', 'Interpreter','latex', FontSize=20);

legend('Meridian lon = 0 deg', 'Equator', 'GT' , 'first point', 'last point', Interpreter='latex', FontSize=20);

end