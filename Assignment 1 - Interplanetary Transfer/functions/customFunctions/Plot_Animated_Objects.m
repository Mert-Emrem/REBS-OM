function [Body, X, Y, Z] = Plot_Animated_Objects(SF, Object, Pos)
% Function to plot and animate celestial objects in 3D space
% This function creates a 3D representation of a celestial body and applies
% a texture map (e.g., a planet's surface) to the plotted sphere.
% 
% PROTOTYPE:
%  [Body, X, Y, Z] = Plot_Animated_Objects(SF, Object, Pos)
% 
% INPUT:
%  SF      [1x1]   Scaling factor for the object's radius [-]
%  Object  [1x1]   Struct containing the celestial object's data:
%                  - Object.Radius: Object's radius [km]
%                  - Object.name:   Name of the object (used for texture)
%  Pos     [3x1]   Position of the object in 3D space [km]
% 
% OUTPUT:
%  Body    [1x1]   Surface handle for the plotted object
%  X       [MxN]   X-coordinates of the sphere's mesh [km]
%  Y       [MxN]   Y-coordinates of the sphere's mesh [km]
%  Z       [MxN]   Z-coordinates of the sphere's mesh [km]
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Scale the object's radius using the provided scaling factor
R = Object.Radius * SF;  

% Generate a sphere with a resolution of 50x50
[X, Y, Z] = sphere(50); 

% Scale the sphere's coordinates to the object's radius
X = X * R; 
Y = Y * R;
Z = Z * R;

% Create a surface object representing the celestial body
Body = surf(Pos(1) + X, Pos(2) + Y, Pos(3) + Z, 'EdgeColor', 'none');

% Read the texture image for the celestial body
image_name = strcat(Object.name, '.jpg'); 
BodyTexture = imread(image_name);

% Apply the texture map to the surface object
Body.CData = flipud(BodyTexture);  % Flip the image to align correctly
Body.FaceColor = 'texturemap';     % Set the texture mapping
Body.FaceAlpha = 1;                % Set full opacity for the object
end
