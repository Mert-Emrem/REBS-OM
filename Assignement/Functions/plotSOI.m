function [] = plotSOI(r_SOI)
% ------------------------------------------------------------------
%
% function [] = plotSOI(r_SOI)
% 
% This function plots the Sphere of Influence of a planet
%
% INPUT:
% r_SOI     [1x1]     [km]       radius of the SOI of a planet
% 
% Author: Serlini Mariagiulia
% Last update: 3/01/2025
%
% ------------------------------------------------------------------

[X_sphere,Y_sphere,Z_sphere] = sphere(50);
X_sphere = X_sphere*r_SOI;
Y_sphere = Y_sphere*r_SOI;
Z_sphere = Z_sphere*r_SOI;
SOI = surf(X_sphere,Y_sphere,Z_sphere);
set(SOI, 'FaceColor', 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'c', 'EdgeAlpha', 0.2); 
hold on

end