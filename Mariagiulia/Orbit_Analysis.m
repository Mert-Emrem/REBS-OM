function [T_ode, Y] = Orbit_Analysis(r0, v0, mu, tspan, type_pb, plot_orb, other_plots, varargin)

% state vector
y0 = [r0; v0]; 

% Set options for the ODE solver
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14); 

if strcmp(type_pb, 'non_perturbed')
        [T_ode, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu), tspan, y0, options);  
        if strcmp(plot_orb, 'no_plot_orbit')
           disp('No orbit plot required')
        elseif strcmp(plot_orb, 'plot_orbit')
            figure
            plot3(Y(:,1), Y(:,2), Y(:, 3), 'k');
            xlabel('X [km]'); 
            ylabel('Y [km]');
            zlabel('Z [km]');
            title('Non perturbed two-body problem orbit');
            axis equal;
            grid on;
        else
            error('Not recognised. Use ''no_plot_orbit'' or ''plot_orbit''.')
        end
    elseif strcmp(type_pb, 'perturbed')
        J2 = varargin{1};
        Re = varargin{2};
        [T_ode, Y] = ode113(@(t,y) ode_2body_pb(t, y, mu, J2, Re), tspan, y0, options);    
        if strcmp(plot_orb, 'no_plot_orbit')
           disp('No orbit plot required')
        elseif strcmp(plot_orb, 'plot_orbit')
            figure
            scatter3(Y(:, 1), Y(:, 2), Y(:, 3), 3, tspan, 'filled');
            colormap('jet')
            colorbar
            xlabel('X [km]'); 
            ylabel('Y [km]');
            zlabel('Z [km]');
            title('Perturbed two-body problem orbit');
            axis equal;
            grid on;
        else
            error('Not recognised. Use ''no_plot_orbit'' or ''plot_orbit''.')
        end
    else
        error('Not recognised. Use ''non_perturbed'' or ''perturbed''.')
end

if strcmp(other_plots, 'no_plots')
           disp('No other plot required')
        elseif strcmp(other_plots, 'more_plots')
            
            % check the constants
            r = Y(:, 1:3)';
            rnorm = vecnorm(r);
            v = Y(:, 4:6)';
            vnorm = vecnorm(v);
            
            % angular momentum
            h_vect_0 = cross(r0, v0);
            h = cross(r, v);
            hnorm = vecnorm(h);
            
            % energy
            en = vnorm.^2./2 - mu./rnorm;
            
            % eccentricity
            e = 1/mu .* cross(v, h) - r./rnorm;
            enorm = vecnorm(e);
            e_h_prod = dot(e, h);
            
            % radial and transversal velocity
            ur = r./rnorm;
            uh = h./hnorm;
            ut = cross(uh, ur);
            
            vr = dot(v, ur);
            vt = dot(v, ut);
            
            
            figure
            plot(T_ode, h(1, :), 'b--', LineWidth=2);
            hold on
            plot(T_ode, h(2, :), 'r--', LineWidth=2);
            plot(T_ode, h(3, :), 'g--', LineWidth=2);
            plot(T_ode, hnorm, 'k--', LineWidth=2);
            grid on
            title('Angular momentum conservation');
            xlabel('Time [s]');
            ylabel('h_x, h_y, h_x, ||h|| [km^2/s]');
            legend('h_x', 'h_y', 'h_z', 'h_{norm}');
            
            figure
            plot(T_ode, en);
            title('Energy conservation');
            xlabel('Time [s]');
            ylabel('\epsilon [km^2/s^2]');
            grid on
            
            figure
            plot(T_ode, e(1, :));
            hold on
            plot(T_ode, e(2, :));
            plot(T_ode, e(3, :));
            plot(T_ode, enorm);
            title('eccentricity')
            legend('e_x', 'e_y', 'e_z', '||e||');
            grid on
              legend('e_x', 'e_y', 'e_z', 'e_{norm}');
            
            figure
            plot(T_ode, rnorm)
            title('position')
            grid on
            
            figure
            plot(T_ode, e_h_prod);
            title('e*h')
            grid on
            
            figure
            plot(T_ode, vr, 'r', LineWidth=2);
            hold on
            plot(T_ode, vt, 'b', LineWidth=2);
            legend('v_r', 'v_{\theta}');
            title('v_r and v_{\theta}');
            grid on
        else
            error('Not recognised. Use ''no_plots'' or ''more_plots''.')
        end



