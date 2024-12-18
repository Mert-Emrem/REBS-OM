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

[~,~,~,~,vt1_i,vf,~,~] = lambertMR(rr_d,rr_f,(t_f-t_d)*24*3600,ksun,0,0,0,0);

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



figure;

hold on
plot3(y1(:,1),y1(:,2),y1(:,3), 'LineWidth', 2)
plot3(y2(:,1),y2(:,2),y2(:,3), 'LineWidth', 2)
plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','- -')
plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','- -')
plot3(r3(1,:), r3(2, :), r3(3, :),'LineStyle','- -')
legend('arc 1', 'arc2','Mercury-Orb', 'Mars-Orb', 'Asteroid N40-Orb')


title('Interplanetary Trajecotry')
grid on

axis equal
hold off

% simulatore(rr_d, vt1_i', vf', (t_f-t_d)*24*3600, ksun, 1)

% x_trail = [];
% y_trail = [];
% z_trail = [];
% 
% % h_point = plot3(y1(1,1),y2(2,1),y1(3,1),'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
% x = y1(1, 1);
% y = y1(1, 2);
% z = y1(1, 3);
% 
% for t = 1:size(y1,2)
% 
%         line = animatedline(x, y, z,'LineWidth',5);
% 
%         x = y1(t, 1);
%         y = y1(t, 2);
%         z = y1(t, 3);
% 
%         addpoints(line, x, y, z);
% 
%         drawnow;
% 
%         % hold on
% 
%         % x_trail = [x_trail, x];
%         % y_trail = [y_trail, y];
%         % z_trail = [z_trail, z];
% 
%         % set(h_point, 'XData', x, 'YData', y, 'ZData', z);
%         % plot3(x_trail, y_trail, z_trail, 'r-'); % Mostra la traccia
% 
%         pause(0.01);
% 
%         if t == size(y1,2)
% 
%             for k = 1:size(y2,2)
% 
%                             x = y2(k, 1);
%                             y = y2(k, 2);
%                             z = y2(k, 3);
% 
%                             line = animatedline(x, y, z,'LineWidth',5);
%                             addpoints(line, x, y, z);
% 
%                             drawnow;
% 
%                             pause(0.01);
%             end
%         end
% 
% end
% plot3(y1(:,1),y1(:,2),y1(:,3))
% plot3(y2(:,1),y2(:,2),y2(:,3))
% plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','-')
% plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','-')
% plot3(r3(1,:), r3(2, :), r3(3, :), 'LineStyle','-')
% legend('Sun','arc 1', 'arc2','Mercury-Orb', 'Mars-Orb', 'Asteroid N40-Orb')


end