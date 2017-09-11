% calculates dimensions of cylindrical RBC scaled to match the experimental
% MCVs. An RBC of average MCV would have reference diameter Lr_ref. Use
% sample number as third parameter to return only the size of that
% particular sample. 
function [Lr, Lz] = calcRBCSize(MCV_, Lr_ref, varargin)

    MCV_ref = mean(MCV_) * 1e-18;
    Lz_ref  = MCV_ref * 4 / pi / Lr_ref / Lr_ref; % reference RBC height
    
    if (nargin > 2)
        
        sample = varargin{1};
        
        % scale RBC so that MCV is matched
        Lr = Lr_ref * (mean(MCV_(sample)) / MCV_ref).^(1/3) * 1e-6;
        Lz = Lz_ref * (mean(MCV_(sample)) / MCV_ref).^(1/3) * 1e-6;
        
    else

        % scale RBC so that MCV is matched
        Lr = Lr_ref * (MCV_ ./ MCV_ref).^(1/3) * 1e-6;
        Lz = Lz_ref * (MCV_ ./ MCV_ref).^(1/3) * 1e-6;
    end

end