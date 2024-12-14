function [Y, Ypert] = OrbitAnalysis(r0, v0, y0, mu_E, J2, Re, tspan)

if nargin<7
% set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
tspan = linspace( 0, 5*T, 1000 );
end

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

dy = @(t,y) ode_2body_pb(t, y, mu_E);
[ T, Y ] = ode113(dy, tspan, y0, options);

if nargin>4
    dy_J2 = @(t,y) ode_2body_pb(t, y, mu_E, J2, Re);
    [T, YPert ] = ode113(dy_J2, tspan, y0, options);
    figure
    scatter3(Y(:, 1), Y(:, 2), Y(:, 3), 0.5, 'r');
    hold on                                                                                          
    scatter3(YPert(:, 1), YPert(:, 2), YPert(:, 3), 3, tspan, 'filled');
    colormap('jet')
    % plot3(YPert(:, 1), YPert(:, 2), YPert(:, 3));
    colorbar
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('Two-body problem orbit');
    axis equal;
    grid on;
    Y = YPert;

else

figure
plot3(Y(:, 1), Y(:, 2), Y(:, 3));
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

end

r = [Y(:, 1), Y(:, 2), Y(:, 3)];
%rnorm = sqrt(r(:, 1).^2 + r(:, 2).^2 + r(:, 3).^2);
rnorm = vecnorm(r, 2, 2);
v = [Y(:, 4), Y(:, 5), Y(:, 6)];

% % angular momentum 
% h = cross(r, v, 2);
% %h_norm = sqrt(h(:, 1).^2 + h(:, 2).^2 + h(:, 3).^2);
% h_norm = vecnorm(h, 2, 2);
% 
% 
% figure
% plot(T, h(:, 1), '--', 'Color','m')
% hold on
% plot(T, h(:, 2), '--', 'Color', 'r')
% plot(T, h(:, 3), '--', 'color' , 'b')
% plot(T, h_norm, 'LineWidth',1.5,'color','k');
% 
% % eccentricity
% e = (1/mu_E).*cross(v, h, 2)- r./rnorm;
% 
% %e_norm = sqrt(e(:, 1).^2 + e(:, 2).^2 + e(:, 3).^2);
% e_norm = vecnorm(e, 2, 2);
% figure
% plot(T, e(:, 1), '--', 'Color','m')
% hold on
% plot(T, e(:, 2), '--', 'Color', 'r')
% plot(T, e(:, 3), '--', 'Color' , 'b')
% %plot(T, e_norm, 'LineWidth',1.5,'Color','k');
% 
% % perpendicular conditon for e and h
% err_edoth =  dot(h, e, 2);
% figure
% plot(T , err_edoth, 'LineWidth',2)
% title('e-h dot product')
% 
% 
% % specific energy
% eps = dot(v, v,2)./2 - mu_E./rnorm;
% figure
% plot(T, eps);
% 
% % radial unit vector
% u_r = r./(rnorm.*ones(size(r)));
% % out of plane vector
% u_h = h./(h_norm.*ones(size(h)));
% % transerval unit vector
% u_t = cross(u_h, u_r);
% % unitary control
% % normut = vecnorm(u_t, 2, 2);
% 
% v_rad = dot(v,u_r, 2);
% v_trv = dot(v,u_t, 2);
% 
% figure
% plot(T, v_rad, 'LineWidth',2)
% hold on
% plot(T, v_trv, 'LineWidth',2)
% 
end