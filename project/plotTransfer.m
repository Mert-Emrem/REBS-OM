function plotTransfer(x, r1, r2, r3)
% 
% Function to plot the two transfer arcs given the times of
% departure, flyby and arrival.
% 
% PROTOTYPE:
%  plotTransfer(x)
% 
% INPUT:
%  x [3,1]        MJD2000 times array: x(1) = departure
%                                      x(2) = flyby
%                                      x(3) = arrival
%  name           Figure title
%  input          Input structure
% 
% OUTPUT:
% plot of three differents arcs

t_d = x(1); t_f = x(2); t_a = x(3); % Initialize times

% Transfer Arcs:

[kep_d,ksun] = uplanet(t_d, 1);
[rr_d,~] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, 4);
[rr_f,~] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,~,~] = ephAsteroids(t_a, 40);
[rr_a,~] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,~,vt1_i,~,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,~,vt2_i,~,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

% Integration of transfer arcs:

[~,y1] = twoBodyInt([0 (t_f-t_d)*24*3600],[rr_d;vt1_i'],ksun);

[~,y2] = twoBodyInt([0 (t_a-t_f)*24*3600],[rr_f;vt2_i'],ksun);


% Plot

figure;

planet = 'Sun';
opts.Units = 'km';
opts.Position = [0, 0, 0];

planet3D(planet, opts);
view([54, 32])
hold on
plot3(y1(:,1),y1(:,2),y1(:,3))
plot3(y2(:,1),y2(:,2),y2(:,3))
plot3(r1(1,:), r1(2, :), r1(3, :))
plot3(r2(1,:), r2(2, :), r2(3, :))
plot3(r3(1,:), r3(2, :), r3(3, :))
legend('Sun','arc 1', 'arc2','Mercury-Orb', 'Mars-Orb', 'Asteroid N40-Orb')


title('Interplanetary Trajecotry')
grid on

axis equal
hold off