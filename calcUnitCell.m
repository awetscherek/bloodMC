% calculates size of unit cell, such that volume fraction of RBC equals HCT
function [L1, L2] = calcUnitCell(Lr, Lz, HCT_)

    L1 = zeros(size(Lr));
    L2 = zeros(size(Lz));
    
    for sample = 1:numel(L1)

        % adjust volume of unit cell accordingly 
        V_cell = @(x) (Lz(sample) + x * 1e-6)*(Lr(sample) + x * 1e-6).^2;
        V_RBC  = Lr(sample).^2 * Lz(sample) * pi / 4;
        L3 = fzero(@(x) (V_cell(x) * HCT_(sample) - V_RBC), 0.4) * 1e-6;
        L1(sample) = Lr(sample) + L3;
        L2(sample) = Lz(sample) + L3;
    end
    
end