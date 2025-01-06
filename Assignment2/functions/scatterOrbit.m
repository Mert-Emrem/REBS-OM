function [] = scatterOrbit(Y, tspan)
% ------------------------------------------------------------------
%
% function [] = scatterOrbit(Y, tspan)
% 
% This function plots the orbitswith scatter given the position vectors at different
% time instants. It requires the tspan compatible with Y as input
%
% INPUTS:
% Y             [nx3]     [km]      matrix of positions vectors in cartesian coordinates: 
%                                   1st col: x components [km]
%                                   2nd col: y components [km]
%                                   3st col: z components [km]
% tspan         [nx1]     [s]       appropriate time span for scatter
%                                   comaptible with Y
% 
% Authors: Serlini Mariagiulia, Bernasconi Ludovico, Emrem Mert, Richero
%         Giovanni
% Last update: 5/01/2025
%
% ------------------------------------------------------------------

scatter3(Y(:,1), Y(:,2), Y(:,3), 2, tspan, 'filled');
colormap('jet')
colorbar
ylabel(colorbar, 'Time [s]', Interpreter='latex', FontSize=14)
hold on

end
