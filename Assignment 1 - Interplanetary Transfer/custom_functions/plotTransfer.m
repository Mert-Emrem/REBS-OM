function plotTransfer(x)
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
% 
% OUTPUT:
% plot of three differents arcs

t_d = x(1); t_f = x(2); t_a = x(3); % Initialize times

% Transfer Arcs:

[kep_d,ksun] = uplanet(t_d, 1);
[rr_d,v_merc] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_f, 4);
[rr_f,v_mars] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,~,~] = ephAsteroids(t_a, 40);
[rr_a,v_harm] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

[~,~,~,~,vt1_i,~,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

[~,~,~,~,vt2_i,~,~,~] = lambertMR(rr_f,rr_a,(t_a-t_f)*24*3600,ksun,0,0,0,0);

% Integration of transfer arcs:

[~,y1] = twoBodyInt([0 (t_f-t_d)*24*3600],[rr_d;vt1_i'],ksun);

[~,y2] = twoBodyInt([0 (t_a-t_f)*24*3600],[rr_f;vt2_i'],ksun);

[~,r1] = twoBodyInt([0 (t_a-t_d)*24*3600],[rr_d;v_merc],ksun);

[~,r2] = twoBodyInt([0 (t_a-t_d)*24*3600],[rr_f;v_mars],ksun);

[~,r3] = twoBodyInt([0 (t_a-t_d+1000)*24*3600],[rr_a;v_harm],ksun);

r1 = r1';
r2 = r2';
r3 = r3';

%% Objects information

Mars.Radius = 3390;
Mars.name = 'Mars';

Mercury.Radius = 2439.7;
Mercury.name = 'Mercury';

Harmonia.Radius = 107.6;
Harmonia.name = 'Harmonia';

Sun.Radius = 700e+3;
Sun.name = 'Sun';

%% Interplanetary transfers plot

figure
hold on
plot3(y1(:,1),y1(:,2),y1(:,3), 'LineWidth', 2, 'Color', '#0714fa')
plot3(y2(:,1),y2(:,2),y2(:,3), 'LineWidth', 2, 'Color', '#e607fa')
plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','- -', 'Color', '#8fe866')
plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','- -', 'Color', '#eb8552')
plot3(r3(1,:), r3(2, :), r3(3, :),'LineStyle','- -', 'Color', '#bfb588')

%Sun
plotObjects(10, Sun, [0, 0, 0])
% Mercury
plotObjects(1800, Mercury, r1(:,1))
% Mars 
plotObjects(2000, Mars, r2(:,1))
% Harmonia
plotObjects(50000, Harmonia, r3(:,1))

legend('arc 1', 'arc 2','Mercury Orbit', 'Mars Orbit', 'Harmonia Orbit')


title('Interplanetary Trajectory')
grid on

axis equal
hold off



end