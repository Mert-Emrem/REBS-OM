function [] = plotGroundTrack(lon, lat, t)
% Function plotGroundTrack plots the groundtrack of a given orbit. It
% should be used after function groundTrack, wich give as output the input
% of this function. 
% Notice: upload the image of the planet in the same folder of this
% function.
%
% STRUCTURE:
%   [] = plotGroundTrack(lon, lat, t)
%
% INPUT:
% - lon       [nx1] longitude
% - lat       [nx1] latitude
% - t         [nx1]   column vector of time in witch to compute the solution
%
% AUTHORS: Bernasconi Ludovico, Emrem Mert, Richero Giovanni, Serlini
%          Mariagiulia
%
% -------------------------------------------------------------------------

% Carica l'immagine
img = imread('MercuryTexture.jpg');

% Ottieni le dimensioni dell'immagine
%[rows, cols, ~] = size(img);

% Crea il plot
figure;
ax = axes;

% Mostra l'immagine come sfondo, ridimensionata per coprire l'intero range degli assi
imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', ax);
set(ax, 'YDir', 'normal');
hold on;
axis on;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

% Imposta i limiti degli assi per coprire l'immagine
xlim([-180 180]);
ylim([-90 90]);

% Tracciare una griglia per i meridiani e paralleli (opzionale)
plot([0 0], ylim, 'k', 'LineWidth', 1); % Meridiano di Greenwich
plot(xlim, [0 0], 'k', 'LineWidth', 1); % Equatore

scatter(lon, lat, 5, t, 'filled');
colorbar
plot(lon(1), lat(1), 'or', LineWidth=2);
plot(lon(end), lat(end), 'ob', LineWidth=2);

end