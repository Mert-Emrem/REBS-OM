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

[~,y1] = twoBodyInt([0 (t_f-t_d-100)*24*3600],[rr_d;vt1_i'],ksun);

[~,y2] = twoBodyInt([0 (t_a-t_f)*24*3600],[rr_f;vt2_i'],ksun);

[~,r1] = twoBodyInt([0 88*24*3600],[rr_d;v_merc],ksun);

[~,r2] = twoBodyInt([0 690*24*3600],[rr_f;v_mars],ksun);

[~,r3] = twoBodyInt([0 1260*24*3600],[rr_a;v_harm],ksun);

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


[kep_d,ksun] = uplanet(t_d, 1);
[rr1_0,v_merc0] = kep2car(kep_d(1),kep_d(2),kep_d(3),kep_d(4),kep_d(5),kep_d(6),ksun);
[kep_f,ksun] = uplanet(t_d, 4);
[rr2_0,v_mars0] = kep2car(kep_f(1),kep_f(2),kep_f(3),kep_f(4),kep_f(5),kep_f(6),ksun);
[kep_a,~,~] = ephAsteroids(t_d, 40);
[rr3_0,v_harm0] = kep2car(kep_a(1),kep_a(2),kep_a(3),kep_a(4),kep_a(5),kep_a(6),ksun);

ToF1 = (t_f-t_d)*24*3600;
ToF2 = (t_a-t_f)*24*3600;
dt = 250000;
tspan = 0:dt:ToF1+ToF2;

figure

hold on
grid on
view(45, 25) % azimuth and elevation
axis equal
xlim([-0.39e+9 +0.41e+9])
ylim([-0.39e+9 +0.41e+9])

% 3 orbits
plot3(r1(1,:), r1(2, :), r1(3, :), 'LineStyle','- -', 'Color', '#8fe866')
plot3(r2(1,:), r2(2, :), r2(3, :), 'LineStyle','- -', 'Color', '#eb8552')
plot3(r3(1,:), r3(2, :), r3(3, :),'LineStyle','- -', 'Color', '#bfb588')

% Lambert arcs
line_y1 = animatedline(rr_d(1), rr_d(2), rr_d(3),'color','#0714fa','LineWidth',2);
line_y2 = animatedline(rr_f(1), rr_f(2), rr_f(3),'color', '#e607fa','LineWidth',2);

% 3 object's orbits
line1 = animatedline(rr1_0(1), rr1_0(2), rr1_0(3),'color', '#8fe866','LineWidth', 1.5);
line2 = animatedline(rr2_0(1), rr2_0(2), rr2_0(3),'color', '#eb8552','LineWidth', 1.5);
line3 = animatedline(rr3_0(1), rr3_0(2), rr3_0(3),'color', '#bfb588' ,'LineWidth',1.5);

% SpaceCraft
SpaceCraft = scatter3(rr_d(1), rr_d(2), rr_d(3), 30, 'k', '<', 'filled');

% Sun
plotObjects(15, Sun, [0, 0, 0])
% Mercury
[ObjMerc, Merc_X, Merc_Y, Merc_Z] = Plot_Animated_Objects(2000, Mercury, rr1_0);
% Mars
[ObjMars, Mars_X, Mars_Y, Mars_Z] = Plot_Animated_Objects(2000, Mars, rr2_0);
% Harmonia
[ObjHarm, Harm_X, Harm_Y, Harm_Z] = Plot_Animated_Objects(50000, Harmonia, rr3_0);


pause(3.5)

r1 = [];
r2 = [];
r3 = [];
y = [];


for jj=1:length(tspan)-1

        if tspan(jj) <= ToF1
        [~,y_temp] = twoBodyInt([0    tspan(jj+1)],[rr_d;vt1_i'],ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y1, y(1,1), y(1,2), y(1, 3));
        else
        [~,y_temp] = twoBodyInt([ToF1 tspan(jj+1)],[rr_f;vt2_i'],ksun);
        y = [y_temp(end, 1:3); y];
        addpoints(line_y2, y(1,1), y(1,2), y(1, 3));
        end

        
        [~,rtemp1] = twoBodyInt([0 tspan(jj+1)],[rr1_0;v_merc0],ksun);
        r1 = [rtemp1(end, 1:3); r1];
        [~,rtemp2] = twoBodyInt([0 tspan(jj+1)],[rr2_0;v_mars0],ksun);
        r2 = [rtemp2(end, 1:3); r2];
        [~,rtemp3] = twoBodyInt([0 tspan(jj+1)],[rr3_0;v_harm0],ksun);
        r3 = [rtemp3(end, 1:3); r3];
        

        set(ObjMerc, 'XData', r1(1,1) + Merc_X, 'YData', r1(1,2) + Merc_Y, 'ZData', r1(1,3) + Merc_Z);
        set(ObjMars, 'XData', r2(1,1) + Mars_X, 'YData', r2(1,2) + Mars_Y, 'ZData', r2(1,3) + Mars_Z);
        set(ObjHarm, 'XData', r3(1,1) + Harm_X, 'YData', r3(1,2) + Harm_Y, 'ZData', r3(1,3) + Harm_Z);
        
        set(SpaceCraft, 'XData', y(1,1), 'YData', y(1,2), 'ZData', y(1, 3))
        
        addpoints(line1, r1(1, 1), r1(1, 2), r1(1, 3));
        addpoints(line2, r2(1, 1), r2(1, 2), r2(1, 3));
        addpoints(line3, r3(1, 1), r3(1, 2), r3(1, 3));
        drawnow;
       
        pause(0.0000001)

end
legend('Mercury Orbit', 'Mars Orbit', 'Harmonia Orbit', 'Lambert arc (1)', 'Lambert arc (2)')


title('Interplanetary Trajectory')

hold off

end

