function Mercury_3D(R_Me)
% -----------------------------------------------------------------------
% Mercury_3D(R_Me)
%
% Plot planet Mercury
% 
%
% INPUT:
%   Rt          [1x1]   [km]       Mercury radius       
%
% OUTPUT:
%   []          [plot]  [-]        Plot Mercury
%
% AUTHORS: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%          Giovanni
% Last Update: 5/1/2025
%
% ------------------------------------------------------------------------

% Set Mercury radius from astroconstant function in case of no inputs 
if nargin < 1
    R_Me = astroConstants(21);                                      
end

% Mercury texture
Mercury_image = 'MercuryTexture.jpg';

% background color: white
background_plot = 'w';

figure('Color', background_plot);
hold on;
grid on;
axis equal;

xlabel('X [km]', Interpreter='latex', FontSize=12);
ylabel('Y [km]', Interpreter='latex', FontSize=12);
zlabel('Z [km]', Interpreter='latex', FontSize=12);

% Initial view
view(120,30);

% Define sphere to create the planet 
npanels = 180;  
[x, y, z] = ellipsoid(0, 0, 0, R_Me, R_Me, R_Me, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

% Load Mercuty image 
cdata = imread(Mercury_image);

% Set transparency = opaque
alpha = 1; 

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end

