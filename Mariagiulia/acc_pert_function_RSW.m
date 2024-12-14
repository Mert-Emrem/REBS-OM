function acc_pert_vect = acc_pert_function_RSW(~, kep_el, J2, mu_E, Re)

[r, ~] = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);

rnorm = norm(r);

i = kep_el(3);
th = kep_el(6);
om = kep_el(5);


a_J2_SRW = -3/2 * J2 * mu_E * Re^2 / rnorm^4 * [1-3*(sin(i))^2*(sin(th+om))^2; 
                                                (sin(i))^2*(sin(2*(th+om))); 
                                                sin(2*i)*sin(th+om)];

acc_pert_vect = a_J2_SRW;

end
