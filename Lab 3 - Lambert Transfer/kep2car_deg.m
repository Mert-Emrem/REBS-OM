function [r_ECI, v_ECI] = kep2car_deg(varargin)
    
% Accepts following sets of arguments:
% - [1] Keplerian vector
% - [2] Keplerian vector, mu
% - [6] Keplerian values
% - [7] Keplerian values, mu

switch nargin
    case 1
        % Single argument: [a, e, i, OM, om, th]
        KEP = varargin{1};
        [a, e, i, OM, om, th] = deal(KEP(1), KEP(2), ...
                                deg2rad(KEP(3)), deg2rad(KEP(4)), ...
                                deg2rad(KEP(5)), deg2rad(KEP(6)));
        mu = astroConstants(13); % Use default gravitational parameter
    case 2
        % Two arguments: [a, e, i, OM, om, th], mu
        KEP = varargin{1};
        [a, e, i, OM, om, th] = deal(KEP(1), KEP(2), ...
                                deg2rad(KEP(3)), deg2rad(KEP(4)), ...
                                deg2rad(KEP(5)), deg2rad(KEP(6)));
        mu = varargin{2}; % Use provided gravitational parameter
    case 6
        % Six arguments: a, e, i, OM, om, th (no mu provided)
        [a, e, i, OM, om, th] = deal(varargin{1:6});
        [i, OM, om, th] = deal(deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th));
        mu = astroConstants(13); % Use default gravitational parameter
    case 7
        % Seven arguments: a, e, i, OM, om, th, mu
        [a, e, i, OM, om, th] = deal(varargin{1:6});
        [i, OM, om, th] = deal(deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th));
        mu = varargin{7};
    otherwise
            error('Invalid number of input arguments. Expecting 1, 2, 6, or 7 inputs.');
end

p = a*(1-e^2);
h = sqrt(p*mu);
r = p / (1+e*cos(th));

r_PF = r*[cos(th), sin(th), 0]';
v_PF = (mu/h) * [-sin(th), (e+cos(th)), 0]';

% Rotation matrices: Earth-Centered Inertial --> Perifocal   (ECI->PF)
R_om = [cos(om)  sin(om)    0   ;
        -sin(om) cos(om)    0   ;
           0        0       1   ];
       
R_i  = [   1        0       0   ;
           0     cos(i)   sin(i);
           0    -sin(i)   cos(i)];
       
R_OM = [cos(OM)  sin(OM)    0   ;
        -sin(OM) cos(OM)    0   ;
           0        0       1   ];
    
R313 = R_om * R_i * R_OM; % ECI --> PF


% PF --> ECI
r_ECI = R313'* r_PF;
v_ECI = R313'* v_PF;



end

