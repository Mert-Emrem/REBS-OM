function dy = eq_motion_cartesian(t, cart_coord, acc_pert_function_cartesian,  mu) %, J2, Re, day)

r = cart_coord(1:3);
v = cart_coord(4:6);
rnorm = norm(r);


%a_p_RSW = acc_pert_function_RSW(t, kep_el); % ) 
a_p = acc_pert_function_cartesian(t, cart_coord);  

dy = [ v(1);
       v(2);
       v(3);
       (-mu/rnorm^3)*r(1) + a_p(1);
       (-mu/rnorm^3)*r(2) + a_p(2);
       (-mu/rnorm^3)*r(3) + a_p(3)];

end