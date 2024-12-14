function plotOrbitDin_Mul_GraphOpts(vet_a, vet_e, vet_i, vet_OM, vet_om, vet_th_i, vet_th_f, ndth, mu, graphopts)

% --------------------------------------------------------------------
% 3D Orbit plotting
%
% plotOrbitDin_Mul_GraphOpts(vet_a, vet_e, vet_i, vet_OM, vet_om, vet_thi, vet_thf, ndth, mu, graphopts)
% 
% The function returns a figure in which the earth and 
% the n orbits under consideration are 
%
% --------------------------------------------------------------------
% Input Values:
% a       [1xn]    semi major axis                                   [km]
% e       [1xn]    eccentricity                                      [-]
% i       [1xn]    inclination                                       [rad]
% OM      [1xn]    RAAN                                              [rad]
% om      [1xn]    pericenter anomaly                                [rad]
% th0     [1xn]    initial true anomaly                              [rad]
% thf     [1xn]    final true anomaly                                [rad]
% ndth    [1x1]    number of elements of the true anomaly vector     [rad]
% mu      [1x1]    gravitational paramater                           [km^3/s^2]
% graphopts -> struct:
%  - graphopts.linewidth       [1xn]       orbit line width          [int]
%  - graphopts.color           [nx3]       orbit color               [int]
%  - graph_opts.linestyle      [1x1]       sorbit line style         [char]
% --------------------------------------------------------------------


%% Construction of the true anomaly matrix
nn = length(vet_a);

MVetTh = zeros(nn,ndth);

for p = 1:1:nn
    vetTh = linspace(vet_th_i(p),vet_th_f(p),ndth);
    MVetTh(p,1:1:end) = vetTh;
end


%% Plotting
M_pos = zeros(3,ndth);

planet = 'Earth Cloudy';
opts.Units = 'km';

planet3D(planet, opts);
light('Position',[1,0,0]);
view(80,15);
grid on
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
hold on;

for p = 1:1:nn

    for j = 1:1:ndth

        [rr, ~] = par2car(vet_a(p), vet_e(p), vet_i(p), vet_OM(p), vet_om(p), MVetTh(p,j), mu);
        M_pos(1,j) = rr(1);    
        M_pos(2,j) = rr(2);
        M_pos(3,j) = rr(3);

    end

    Line = animatedline('Color', [graphopts.color(p,1) graphopts.color(p,2) graphopts.color(p,3)],...
        'LineWidth', graphopts.linewidth(p), 'LineStyle', graphopts.linestyle);

    for t = 1:1:ndth

        x = M_pos(1,t);
        y = M_pos(2,t);
        z = M_pos(3,t);

        addpoints(Line, x, y, z);

        drawnow;

        pause(0.1);

    end

end


