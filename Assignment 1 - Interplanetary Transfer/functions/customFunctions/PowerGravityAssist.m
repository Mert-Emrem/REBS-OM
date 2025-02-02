function [dvp, t_SOI, rp, delta, e_m, e_p, am, ap, vpm, vpp, deltam, deltap] = PowerGravityAssist(vinfm, vinfp, Rpl, hatm, k_f, dynamic_ON)

% Brief explanation:
% The function calculates the parameters associated with a powered gravity assist, 
% including the turn angles, eccentricities, and velocities for the incoming 
% and outgoing hyperbolic trajectories. The method uses an iterative approach 
% to determine the pericenter radius that satisfies the required turn angle.
%
% 
% PROTOTYPE:
%  [dvp, delta, rp, em, ep, am, ap, vpm, vpp, deltam, deltap] = PowerGravityAssist(vinfm, vinfp, Rpl, hatm, k_f, dynamic_ON)
%  
% INPUT:
%  vinfm [3]       Incoming relative velocity to planet vector [km/s]
%  vinfp [3]       Outgoing relative velocity to planet vector [km/s]
%  Rpl [1]         Radius of the planet [km]
%  hatm [1]        Height of the planet's atmosphere [km]
%  k_f [1]         Gravitational parameter of the planet [km^3/s^2]
%  dynamic_ON [1]  Boolean flag to enable dynamic bounds adjustment
% 
% OUTPUT:
%  dvp [1]         Delta-V required at pericenter [km/s] 
%  delta [1]       Turn angle [rad]
%  rp [1]          Radius of pericenter [km]
%  e_m [1]          Eccentricity of the incoming hyperbola 
%  e_p [1]          Eccentricity of the outgoing hyperbola
%  am [1]          Semi-major axis of the incoming hyperbola [km]
%  ap [1]          Semi-major axis of the outgoing hyperbola [km]
%  vpm [1]         Velocity at pericenter of the incoming hyperbola [km/s]
%  vpp [1]         Velocity at pericenter of the outgoing hyperbola [km/s]
%  deltam [1]      Turn angle of the incoming hyperbola [rad]
%  deltap [1]      Turn angle of the outgoing hyperbola [rad]
%  t_SOI [1]       Time spent inside SOI [s]
% 
% AUTHORS: Richero Giovanni, Emrem Mert, Bernarsconi Ludovico, Serlini Mariagiulia
% -------------------------------------------------------------------------

    % Turn angle between incoming and outgoing hyperbolas
    delta = acos(dot(vinfm, vinfp) / norm(vinfm) / norm(vinfp));
    
    % Eccentricity and turn angle as functions of rp
    emp = @(rp, vinf) 1 + rp * vinf^2 / k_f;
    deltamp = @(rp, vinf) 2 * asin(1 / emp(rp, vinf));
    

    % Nonlinear equation to solve
    f = @(rp) delta - 0.5 * deltamp(rp, norm(vinfm)) - 0.5 * deltamp(rp, norm(vinfp));
    
    % SOI radius
    a_mars = 227.956e6; % [km]
    mu_S = astroConstants(4);

    r_SOI = a_mars*(k_f/mu_S)^(2/5);

    % Minimum radius (planet radius + atmosphere height)
    rp_min = Rpl + hatm;
    upper_bound = r_SOI;

    % Check initial values at the boundaries
    f_min = f(rp_min);
    f_max = f(upper_bound);
    
    rp = rp_min; 

    if dynamic_ON
        % Dynamic interval adjustment
        while isnan(f_max) || ~isreal(f_max) || (f_min * f_max > 0)
            upper_bound = upper_bound / 10;
            f_max = f(upper_bound);
            if upper_bound < 10 * rp_min  % Stop if bounds shrink too much
                break;
            end
        end
        
        % Try to find the solution using fzero
        options = optimset('TolX', 1e-14, 'FunValCheck', 'on');
        try
            rp = fzero(f, [rp_min, upper_bound], options);
        catch
            rp = rp_min;  % Penalize by setting rp to minimum possible value
        end

    elseif (f_min * f_max < 0)
        options = optimset('TolX', 1e-14);
        rp = fzero(f, [rp_min, upper_bound], options);
    end
    
    % Calculate velocities at pericenter
    vpp = sqrt(norm(vinfp)^2 + 2 * k_f / rp);
    vpm = sqrt(norm(vinfm)^2 + 2 * k_f / rp);
    dvp = norm(vpp - vpm);
    
    % Eccentricities of the hyperbolae
    e_m = emp(rp, norm(vinfm));
    e_p = emp(rp, norm(vinfp));
    
    % Turn angles
    deltam = deltamp(rp, norm(vinfm));
    deltap = deltamp(rp, norm(vinfp));
    
    % Semi-major axes of incoming/outgoing hyperbolas
    am = abs(rp / (1 - e_m));
    ap = abs(rp / (1 - e_p));

    % Angular momenta of the hyperbolae
    h_m = vpm * rp;
    h_p = vpp * rp;
    th_m = acos((h_m^2/(k_f*r_SOI)-1)/e_m);
    th_p = acos((h_p^2/(k_f*r_SOI)-1)/e_p);
    
    % Hyperbolic anomalies
    F_th_m = 2 * atanh(sqrt((e_m - 1) / (e_m + 1)) * tan(th_m / 2));
    F_th_p = 2 * atanh(sqrt((e_p - 1) / (e_p + 1)) * tan(th_p / 2));
    
    % Mean anomalies
    M_th_m = e_m * sinh(F_th_m) - F_th_m;
    M_th_p = e_p * sinh(F_th_p) - F_th_p;
    
    % Time inside SOI
    t_SOI_m = (h_m^3 / k_f^2) * M_th_m / (e_m^2 - 1)^(3/2);
    t_SOI_p = (h_p^3 / k_f^2) * M_th_p / (e_p^2 - 1)^(3/2);
    
    t_SOI = t_SOI_m + t_SOI_p;

end
