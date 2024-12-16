function plot_grounfTrack(lat,lon,flag,t0,tf)
% 
% This function plot the ground Track of a satellite.


switch flag
    case 1
        % Upload image
        img = imread('EarthTexture.jpg');

        % Mostra l'immagine come sfondo, ridimensionata per coprire l'intero range degli assi
        imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', gca);
        set(gca, 'YDir', 'normal'); %Invert y axis
        hold on

        % Imposta i limiti degli assi per coprire l'immagine
        xlim([-180 180]);
        ylim([-90 90]);

        % Tracciare una griglia per i meridiani e paralleli (opzionale)
        plot([0 0], ylim, 'r--', 'LineWidth', 2); % Meridiano di Greenwich
        plot(xlim, [0 0], 'r--', 'LineWidth', 2); % Equatore
        plot(lon,lat,'y', 'LineStyle','none','Marker','.')
        axis on;

        % Define initial and final point and plot them
        x_i     = lon(1);
        y_i     = lat(1);
        x_f     = lon(end);
        y_f     = lat(end);

        plot(x_i, y_i, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        plot(x_f, y_f, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r');


        % Aggiungi il nome del punto
        text(x_i, y_i, 'Initial point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        text(x_f, y_f, 'Final point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

        % Aggiungi etichette e titolo
        xlabel('Longitutde');
        ylabel('Latitude');
        title('Ground Track');
        
    case 2
        % Upload image
        img = imread('EarthTexture.jpg');

        % Mostra l'immagine come sfondo, ridimensionata per coprire l'intero range degli assi
        imshow(img, 'XData', [-180 180], 'YData', [90 -90], 'Parent', gca);
        set(gca, 'YDir', 'normal'); %Invert y axis
        hold on

        % Imposta i limiti degli assi per coprire l'immagine
        xlim([-180 180]);
        ylim([-90 90]);

        % Tracciare una griglia per i meridiani e paralleli (opzionale)
        plot([0 0], ylim, 'r--', 'LineWidth', 2); % Meridiano di Greenwich
        plot(xlim, [0 0], 'r--', 'LineWidth', 2); % Equatore
        scatter(lon, lat ,2,linspace(t0,tf,length(lon))); 
        colormap("default")
        colorbar; 
        axis on;

        % Define initial and final point and plot them
        x_i     = lon(1);
        y_i     = lat(1);
        x_f     = lon(end);
        y_f     = lat(end);

        plot(x_i, y_i, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        plot(x_f, y_f, 'r*', 'MarkerSize', 10, 'MarkerFaceColor', 'r');


        % Aggiungi il nome del punto
        text(x_i, y_i, 'Initial point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        text(x_f, y_f, 'Final point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

        % Aggiungi etichette e titolo
        xlabel('Longitutde');
        ylabel('Latitude');
        title('Ground Track');

end


