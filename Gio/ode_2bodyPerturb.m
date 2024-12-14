function dy = ode_2bodyPerturb( t, y, mu, method )

% J2 Constants
    R_e = 6378.137;
    J2 = astroConstants(9);



switch method
    case 'cart'
        rr = y(1:3);
        vv = y(4:6);
        
    case 'gauss'
        a  = y(1);   e  = y(2);   i  = y(3);
        OM = y(4);   om = y(5);   th = y(6);

        [rr, vv] = kep2car(a,e,i,OM,om,th, mu);
        
end

r = norm(rr);
v = norm(vv);
x = rr(1);
y = rr(2);
z = rr(3);


% J2 perturbing acceleration in cartesian coordinates:

    kJ2 = 1.5*J2*mu*R_e^2/r^4;
    a_J2_x = kJ2 * x/r*(5*z^2/r^2-1);
    a_J2_y = kJ2 * y/r*(5*z^2/r^2-1);
    a_J2_z = kJ2 * z/r*(5*z^2/r^2-3);



% acc perturbing ... cartesian coordinates:
    a_p_x = a_J2_x ;
    a_p_y = a_J2_y ;
    a_p_z = a_J2_z ;


if strcmp(method, 'cart')

    dy = [  vv(1)                   ;
            vv(2)                   ;
            vv(3)                   ;
            -mu/r^3 * x  +  a_p_x  ;
            -mu/r^3 * y  +  a_p_y  ;
            -mu/r^3 * z  +  a_p_z  ];

elseif strcmp(method, 'gauss')
    % J2 {t,n,h} coordinates:
    tt = vv/norm(vv);                   % Tangent  unit vector
    hh = cross(rr,vv); 
    hh = hh/norm(hh);                   % Normal   unit vector
    nn = cross(hh,tt);                  % Binormal unit vector
    ROT_tnh2xyz = [tt(:) nn(:) hh(:)];
    a_p_xyz = [a_p_x a_p_y a_p_z]';
     
    a_p_tnh = ROT_tnh2xyz' * a_p_xyz;
    a_t = a_p_tnh(1);  a_n = a_p_tnh(2);  a_h = a_p_tnh(3);
    

    b = a * sqrt(1-e^2);
    p = b^2/a;
    n = sqrt(mu/a^3);
    h = n*a*b;
    th_star = th + om;

    a_dot  = 2*a^2*v/mu * a_t;
    e_dot  = 1/v* ( 2*(e+cos(th))*a_t - r/a*sin(th)*a_n );
    i_dot  = r*cos(th_star)/h * a_h;
    OM_dot = r*sin(th_star)/(h*sin(i)) * a_h;
    om_dot = 1/(e*v) * ( 2*sin(th)*a_t + (2*e + r/a*cos(th))*a_n ) - r*sin(th_star)*cos(i)/(h*sin(i))*a_h;
    th_dot = h/(r^2) - 1/(e*v) * ( 2*sin(th)*a_t + (2*e + r/a*cos(th))*a_n );

    dy = [  a_dot ;
            e_dot ;
            i_dot ;
            OM_dot;
            om_dot;
            th_dot];
end

    
  
 
end
