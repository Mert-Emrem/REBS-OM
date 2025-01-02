function plotObjects(SF, Object, Pos)
% Function to plot celestial objects in 3D space
% This function generates a 3D representation of a celestial object and 
% applies a texture map to create a realistic visualization.
% 
% PROTOTYPE:
%  plotObjects(SF, Object, Pos)
% 
% INPUT:
%  SF      [1x1]   Scaling factor for the object's radius [-]
%  Object  [1x1]   Struct containing the celestial object's data:
%                  - Object.Radius: Object's radius [km]
%                  - Object.name:   Name of the object (used for texture)
%  Pos     [3x1]   Position of the object in 3D space [km]
% 
% OUTPUT:
%  A rendered plot of the celestial object in the current figure window.
%
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

% Scale the object's radius using the provided scaling factor
R = Object.Radius * SF;  

% Generate a sphere with a resolution of 50x50
[X, Y, Z] = sphere(50); 

% Scale the sphere's coordinates to match the object's radius
X = X * R; 
Y = Y * R;
Z = Z * R;

% Create a surface object for the celestial body at the specified position
Body = surf(Pos(1) + X, Pos(2) + Y, Pos(3) + Z, 'EdgeColor', 'none');

% Load the texture image for the celestial object
image_name = strcat(Object.name, '.jpg'); 
BodyTexture = imread(image_name);

% Apply the texture map to the surface object
Body.CData = flipud(BodyTexture);  % Flip the image for correct orientation
Body.FaceColor = 'texturemap';     % Set the texture mapping
Body.FaceAlpha = 1;                % Set full opacity
end
