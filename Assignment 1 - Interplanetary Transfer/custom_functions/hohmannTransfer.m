mu = astroConstants(4);

a1 = 57.909e6; % Mercury
a2 = 227.956e6; % Mars

H1 = hohmannTransferDeltaV(a1, a2, mu);

[kep, f, ~] = ephAsteroids_vec(arr_window, 40);

a3 = kep(1);

H2 = hohmannTransferDeltaV(a2, a3, mu);

Hohmann = H1 + H2


function deltaV = hohmannTransferDeltaV(a1, a2, mu)

    v1 = sqrt(mu / a1);

    v2 = sqrt(mu / a2);

    v_transfer1 = sqrt(2 * mu * a2 / (a1 * (a1 + a2)));

    v_transfer2 = sqrt(2 * mu * a1 / (a2 * (a1 + a2)));

    deltaV1 = abs(v_transfer1 - v1); 
    deltaV2 = abs(v2 - v_transfer2); 

    deltaV = deltaV1 + deltaV2;
end
