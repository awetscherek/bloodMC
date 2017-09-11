% calculates number of water molecules that start inside RBC (rest starts
% in plasma).
function N_e = calcNinRBC(NRep, c_e, c_p, HCT_)

    f_e = c_e * HCT_ ./ (c_e * HCT_ + c_p * (1 - HCT_));

    N_e = round(NRep * f_e);

end