function [] = plotGroundTrack(lon, lat, t)

% Carica l'immagine
img = imread('EarthTexture.jpg');

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