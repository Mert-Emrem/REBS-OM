function Terra
Rt = 6371.01;
Earth_image = "Earth_image.jpeg";
% Choose the color of the figure background
background_plot = 'w';
%x = randn(100, 1);
%y = randn(100, 1);
%z = randn(100, 1);


% Create the figure
figure %('Color', [0.0706 0.1608 0.3020]);
hold on

% Modifica del colore degli assi
ax = gca; % Ottieni il handle dell'oggetto 'axes'
%ax.XColor = 'white'; % Colore dell'asse x
%ax.YColor = 'white'; % Colore dell'asse y
%ax.ZColor = 'white'; % Colore dell'asse z

% Modifica del colore della griglia
grid on; % Attiva la griglia
ax.GridColor = [0.5, 0.5, 0.5]; % Imposta il colore della griglia a grigio chiaro


%scatter3(x, y, z, 'filled');
%set(gca, 'Color', 'k');
%set(gca, 'XColor', 'w', 'YColor', 'w','ZColor', 'w');


% Set the axes scale equal
axis equal;

% Put the axes labels
xlabel('X [km]')% 'Color', 'w');
ylabel('Y [km]')% 'Color', 'w');
zlabel('Z [km]')% 'Color', 'w');

% Set initial view
view(120,30);
% Define the number of panels to be used to model the sphere 
npanels = 180;  


% Create the globe with the surf function
cdata = imread(Earth_image);

[X, Y, Z] = sphere(100);
X = X*Rt;
Y = Y*Rt;
Z = Z*Rt;
surf(X, Y, -Z, 'FaceColor', 'texturemap', 'CData', cdata, 'EdgeColor', 'none');
