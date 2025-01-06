function [] = plotOrbit(Y, line_width)
% ------------------------------------------------------------------
%
% function [] = plotOrbit(Y, line_width)
% 
% This function plots the orbits given the position vectors at different
% time instants and an appropriate line width
%
% INPUTS:
% Y             [nx3]     [km]      matrix of positions vectors in cartesian coordinates: 
%                                   1st col: x components [km]
%                                   2nd col: y components [km]
%                                   3st col: z components [km]
% line_width    [1x1]     [-]       appropriate line width for the plot
% 
% Authors: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 3/01/2025
%
% ------------------------------------------------------------------

plot3(Y(:,1), Y(:,2), Y(:,3), LineWidth=line_width);
hold on

end
