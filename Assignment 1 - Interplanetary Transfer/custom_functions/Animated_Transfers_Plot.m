function Animated_Transfers_Plot(x)
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

[~,r3] = twoBodyInt([0 (t_a-t_d)*24*3600],[rr_a;v_harm],ksun);

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


% figure
% hold on
% plot3(y1(:,1),y1(:,2),y1(:,3), 'LineWidth', 2, 'Color', '#0714fa')
% plot3(y2(:,1),y2(:,2),y2(:,3), 'LineWidth', 2, 'Color', '#e607fa')
% plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','- -', 'Color', '#8fe866')
% plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','- -', 'Color', '#eb8552')
% plot3(r3(1,:), r3(2, :), r3(3, :),'LineStyle','- -', 'Color', '#bfb588')
% 
% %Sun
% plotObjects(10, Sun, [0, 0, 0])
% % Mercury
% plotObjects(1800, Mercury, r1(:,1))
% % Mars 
% plotObjects(2000, Mars, r2(:,1))
% % Harmonia
% plotObjects(50000, Harmonia, r3(:,1))
% 
% legend('arc 1', 'arc 2','Mercury Orbit', 'Mars Orbit', 'Harmonia Orbit')
% 
% 
% title('Interplanetary Trajectory')
% grid on
% 
% axis equal
% hold off

%%

[kep_d,ksun] = uplanet(t_d, 1);
[rr1_0,v_merc0] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_d, 4);
[rr2_0,v_mars0] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,~,~] = ephAsteroids(t_d, 40);
[rr3_0,v_harm0] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

ToF1 = (t_f-t_d)*24*3600;
ToF2 = (t_a-t_f)*24*3600;
dt = 100000;
tspan = [0:dt:ToF1, ToF1+1:dt:ToF2, ToF2];

figure
hold on
plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','- -', 'Color', '#8fe866')
plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','- -', 'Color', '#eb8552')
plot3(r3(1,:), r3(2, :), r3(3, :),'LineStyle','- -', 'Color', '#bfb588')

%Sun
plotObjects(10, Sun, [0, 0, 0])

[ObjMerc, Merc_X, Merc_Y, Merc_Z] = Plot_Animated_Objects(1800, Mercury, rr1_0);
[ObjMars, Mars_X, Mars_Y, Mars_Z] = Plot_Animated_Objects(2000, Mars, rr2_0);
[ObjHarm, Harm_X, Harm_Y, Harm_Z] = Plot_Animated_Objects(50000, Harmonia, rr3_0);


r1 = [];
r2 = [];
r3 = [];
y = [];

line1=animatedline(rr1_0(1), rr1_0(2), rr1_0(3),'color', 'r','LineWidth',2);
line2=animatedline(rr2_0(1), rr2_0(2), rr2_0(3),'color', 'r','LineWidth',2);
line3=animatedline(rr3_0(1), rr3_0(2), rr3_0(3),'color', 'r','LineWidth',2);
line_y=animatedline(rr_d(1), rr_d(2), rr_d(3),'color', 'r','LineWidth',2);

for jj=1:length(tspan)-1
        if tspan<= ToF1
        [~,y_temp] = twoBodyInt([0    tspan(jj+1)],[rr_d;vt1_i'],ksun);
        else
        [~,y_temp] = twoBodyInt([ToF1 tspan(jj+1)],[rr_f;vt2_i'],ksun);
        end
        y = [y_temp(end, 1:3); y];
        
        [~,rtemp1] = twoBodyInt([0 tspan(jj+1)],[rr1_0;v_merc0],ksun);
        r1 = [rtemp1(end, 1:3); r1];
        [~,rtemp2] = twoBodyInt([0 tspan(jj+1)],[rr2_0;v_mars0],ksun);
        r2 = [rtemp2(end, 1:3); r2];
        [~,rtemp3] = twoBodyInt([0 tspan(jj+1)],[rr3_0;v_harm0],ksun);
        r3 = [rtemp3(end, 1:3); r3];
        set(ObjMerc, 'XData', r1(1,1) + Merc_X, 'YData', r1(1,2) + Merc_Y, 'ZData', r1(1,3) + Merc_Z);
        set(ObjMars, 'XData', r2(1,1) + Mars_X, 'YData', r2(1,2) + Mars_Y, 'ZData', r2(1,3) + Mars_Z);
        set(ObjHarm, 'XData', r3(1,1) + Harm_X, 'YData', r3(1,2) + Harm_Y, 'ZData', r3(1,3) + Harm_Z);

        addpoints(line_y, y(1,1), y(1,2), y(1, 3));
        addpoints(line1, r1(1, 1), r1(1, 2), r1(1, 3));
        addpoints(line2, r2(1, 1), r2(1, 2), r2(1, 3));
        addpoints(line3, r3(1, 1), r3(1, 2), r3(1, 3));
        drawnow;
       
        pause(0.00001)

end
axis equal

% 
%         BodyTexture = imread('EarthTexture.jpg');
%         Body.CData = flipud(BodyTexture); % Flip to align texture correctly
%         Body.FaceColor = 'texturemap';
%         Body.FaceAlpha = 1; 
% 
%         hold on
%         plot3(0, 0, 0)
%         pause(0.01)
% 
%         [~,r1] = twoBodyInt([dt(jj), dt(jj+1)],[r1(e);v_merc],ksun);
% 
% 
%             tVect=[t0:delta:t0+dt(jj), t0+dt(jj)];
%                 for kk=1:length(tVect)-1
%                     [t, y]=ode113(g, tVect(kk:kk+1), y0, opts);
%                     addpoints(line, y(:,4), y(:,5), y(:,6));
%                     drawnow;
%                     pause(dtime);
%                     y0=y(end,:);
% 
%                 end
% 
% 
% for t = 1:1:size()
% 
%     % Update Earth Position
% 
% 
%     % Plot Center Point (Optional)
%     plot3(0, 0, 0, 'ro');
% 
%     pause(0.01)

end

