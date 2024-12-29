function plotObjects(SF, Object, Pos)

R = Object.Radius*SF;  % Object radius in Km
        [X, Y, Z] = sphere(50); % 50x50 resolution
        X = X * R; % Scale to object's radius
        Y = Y * R;
        Z = Z * R;

Body = surf(Pos(1)+X, Pos(2)+Y, Pos(3)+Z, 'EdgeColor', 'none');


image_name = strcat(Object.name, '.jpg');

BodyTexture = imread(image_name);
Body.CData = flipud(BodyTexture); % Flip to align texture correctly
Body.FaceColor = 'texturemap';
Body.FaceAlpha = 1; 

