function acc_pert_vect = acc_pert_function_RSW(~, kep_el, J2, mu_E, Re)

[r, ~] = par2car(kep_el(1), kep_el(2), kep_el(3), kep_el(4), kep_el(5), kep_el(6), mu_E);

rnorm = norm(r);

a_J2_SRW = -3/2 * J2 * mu_E * Re / rnorm^4 * [1-3*(sin(kep_el(3)))^2*(sin(kep_el(6)-kep_el(5)))^2; (sin(kep_el(3)))^2*(sin(2*(kep_el(6)-kep_el(5)))); (sin(kep_el(3)))^2*(sin(2*(kep_el(6)-kep_el(5))))];

acc_pert_vect = a_J2_SRW;

end
