% calculates the intervals needed to cover the whole range of Ts (as
% specified in T_MP_) for both gradient profiles.
function [NTimeIntervals, TInterval] = calcIntervalDurations(TE_, T_MP_, varargin)

% last parameter allows to specify a T_FC different from T_MP (as in
% TE-dependence).
if (nargin > 2)
    T_FC_ = varargin{1};
    if (nargin > 3)
        TE_min_FC = varargin{2};
    else
        TE_min_FC = 0;
    end
else
    T_FC_ = T_MP_;
    TE_min_FC = 0;
end

TTT = 0;

% intervals needed for bipolar gradients:
for T = T_MP_
  for TE = TE_

    % The time intervals needed for bipolar diffusion time T_MP:
    TT_ = [(max(TE_) - TE) / 2, (max(TE_) - T) / 2, ...
           [1/2, 1] * T + (max(TE_) - T) / 2, (max(TE_) + TE) / 2]; 

    for TT = TT_
      if (min(abs(TTT - TT)) > 1e-10)
        TTT = [TTT; TT]; %#ok<AGROW>
      end
    end
  end
end

% intervals needed for flow-comp gradients:
for T = T_FC_
  for TE = TE_(TE_ > TE_min_FC)
    % The additional time intervals for flow-comp diffusion time T_FC:
    TT_ = [(max(TE_) - T) / 2, [(2-sqrt(2))/4, 1/2, (2+sqrt(2))/4, 1] ...
           * T + (max(TE_) - T) / 2];
    for TT = TT_
      if (min(abs(TTT - TT)) > 1e-10)
        TTT = [TTT; TT]; %#ok<AGROW>
      end
    end
  end
end

TTT = sort(TTT(2:numel(TTT)));

TInterval = TTT - [0; TTT(1:(numel(TTT)-1))];
NTimeIntervals = numel(TInterval);

end
