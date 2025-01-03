function [] = scatterOrbit(Y, tspan)

scatter3(Y(:,1), Y(:,2), Y(:,3), 2, tspan, 'filled');
colormap('jet')
colorbar
ylabel(colorbar, 'Time [s]')
hold on

end
