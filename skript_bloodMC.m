
% close pool before mex compilation
if ~isempty(gcp('nocreate'))
     delete(gcp('nocreate'));
end
mex MCKernel3D_inout.cpp

%% simulation
clc
clear

fprintf('\nMonte-Carlo Simulation for water diffusion in blood');
fprintf('\n--------------------------------------------------------\n');

settings.NRep    =    300000; % number of simulated particles in manuscript
settings.NRep    =      1000; % for a quick & dirty test that doesn't run
                              % for several days.

settings.D_p     =   2.75e-9; % diffusion in plasma
settings.D_e     =   1.00e-9; % diffusion in RBC (< diffusion in plasma)

settings.ti_ex   =  12.00e-3; % intra-cellular pre-exchange lifetime

settings.c_p     =   0.95;    % free water concentration in Plasma
settings.c_e     =   0.70;    % free water concentration in RBC

settings.Ld_ref  =   8.00e-6; % reference RBC diameter

settings.dr_rel  =   1.00e-2; % desired max step size for MC calculations
                              % relative to shortest dimension, which in 
                              % this case is NOT the RBC height (Lh), but 
                              % the min distance between two neighbouring 
                              % RBCs (L2 - Lh).
  
fprintf('\nSimulation parameters:\n');
disp(settings)
                                                 
% experimental blood sample measurements:
experiment.HCT_ = [37.8 37.8  38.1 39.0 40.9 42.2 43.1 43.7 44.5 47.1]; %%
experiment.MCV_ = [90.6 86.3  83.7 88.6 87.2 83.7 87.4 94.0 85.2 89.9]; %fL
experiment.T2_  = [184. 57.2 125.6 41.7 83.7 37.3 39.6 83.8 77.9 31.8]; %ms

% the editors wanted a more simple figure, so we didn't use the
% experimental values, but only covered the HCT range and set the MCV to
% the mean value
experiment.HCT_ = [37   39   41   43   45   47];
experiment.MCV_ = [87.7 87.7 87.7 87.7 87.7 87.7];

% experimental MRI settings for T-dependence:
experiment.T_        = [40, 50, 60, 70,  80,  90, 100];           %ms
experiment.TE        = 120;                                       %ms

% experimental MRI settings for TE-dependence:
experiment.TE_       = [60, 70, 80, 90, 100, 120, 140, 160, 200]; %ms
experiment.T_MP      =  30;                                       %ms
experiment.T_FC      =  70;                                       %ms

% minimum possible TE for flow-comp with T_FC = 70 ms
experiment.TE_min_FC =  90;                                       %ms

%% water reference:
% DO NOT execute this for blood diffusion simulation!
% overwrites some settings, if you want to do a reference simulation for
% water with extremely high permeable membranes (= free diffusion):

settings.D_p     =   3.0e-9;
settings.D_e     =   3.0e-9;

settings.ti_ex   =    1e-10;  % intra-cellular pre-exchange lifetime

settings.c_e     =      1.0;  % free water concentration in RBC
settings.c_p     =      1.0;  % free water concentration in Plasma

fprintf('\nSimulation parameters:\n');
disp(settings)

experiment.HCT_ = mean(experiment.HCT_);
experiment.MCV_ = mean(experiment.MCV_);
experiment.T2_  = mean(experiment.T2_);

%%
% do the simulation - can take A LOT of time

% open pool before tic (if no pool is there already)
if isempty(gcp('nocreate'))
    [~] = parpool(12);
end

for sample = 1:numel(experiment.HCT_)
    
  fprintf('\n--------------------------------------------------------\n');
    
  [Ld, Lh] = calcRBCSize(experiment.MCV_, settings.Ld_ref, sample);
  [L1, L2] = calcUnitCell(Ld, Lh, experiment.HCT_(sample) * 1e-2);
  
  fprintf(['\nSample No.%2d, geometric HCT = % 3.3f %%, ' ...
      'MCV = % 3.1f fL\n'], sample, 100 * (Ld.^2 .* Lh * pi / 4) ...
      ./ (L1.^2 .* L2), pi * Ld.^2 .* Lh / 4e-18);

  N_e = calcNinRBC(settings.NRep, settings.c_e, settings.c_p, ...
                   experiment.HCT_(sample) * 1e-2);

  fprintf('Particles starting in RBC:    %6d\n', N_e);
  fprintf('Particles starting in Plasma: %6d\n', settings.NRep - N_e);
  
  % permeability = surface-to-volume ratio divided by pre-exchange lifetime
  kappa = Ld * Lh / (2 * Ld + 4 * Lh) / settings.ti_ex;
  
  fprintf(['Membrane permeability for RBC -> plasma ' ...
      '= %2.6f * 10^-3 cm/s\n'], kappa * 1e5);
  
  % Precalculate ratio of kappa and D:
  k2D_p = 2 * kappa / settings.D_p * (settings.c_e / settings.c_p);
  k2D_e = 2 * kappa / settings.D_e;
  
  dr = settings.dr_rel * min(L2 - Lh, Lh);
  
  fprintf('Maximum step size dr = %3.3f nm\n', dr * 1e9);
 
  % calculate time intervals for T-dependence
  [NTimeIntervals, TInterval] = ...
      calcIntervalDurations(experiment.TE * 1e-3, experiment.T_ * 1e-3);
  
  % alternatively: calculate time intervals for TE-dependence
  %[NTimeIntervals, TInterval] = calcIntervalDurations(...
  %    experiment.TE_ * 1e-3, experiment.T_MP * 1e-3, ...
  %    experiment.T_FC * 1e-3, experiment.TE_min_FC * 1e-3);
    
  % minimum number of steps such that the desired maximum step size is
  % not larger than specified in settings.dr
  NSteps = ceil(6 * settings.D_p * TInterval / dr.^2);
  
  fprintf('Total number of steps for each particle = %d\n', sum(NSteps));
        
  % Preparing output matrix
  data = zeros(settings.NRep, NTimeIntervals * 9 + 3);
        
  % Precalculate 2 * step size:
  dr_e = sqrt(24 * settings.D_e * TInterval ./ NSteps);    
  dr_p = sqrt(24 * settings.D_p * TInterval ./ NSteps);

  tic
  % split into 100 parfor loops, so that progress cam be seen:
  for p = 1:100
    
    parfor q = (floor((p-1) * size(data, 1) / 100) + 1) ...
              : min(floor(p * size(data, 1) / 100), size(data, 1))
    
      % saves information about residence times and mean positions:
      tmp = zeros(1, NTimeIntervals * 9 + 3);
            
      % starting position:
      condition = true;
      while (condition)
        tmp(1) = rand(1) * L1;
        tmp(2) = rand(1) * L1;
        tmp(3) = rand(1) * L2;
        condition = xor(((tmp(1)-L1/2)^2 + (tmp(2)-L1/2)^2 > Ld^2/4) ...
                       || ((tmp(3)-L2/2)^2 > Lh^2/4), q > N_e);
      end
        
      for j = 1:NTimeIntervals
          
        % each interval uses 9 fields
        tmp((j*9) - 6 + (1:9)) = MCKernel3D_inout( ...
          tmp(j*9 - 8), tmp(j*9 - 7), tmp(j*9 - 6), L1, L2, Ld/2, Lh/2, ...
          dr_e(j), dr_p(j), k2D_e, k2D_p, NSteps(j)); %#ok<PFBNS>
      end
    
      data(q, :) = tmp;
    end
      
    fprintf('.');
    if (mod(p, 10) == 0)
        fprintf('\n');
    end
  end  
  toc;
    
  if (settings.D_p < 2.9e-9)
    save(sprintf('s%02d_T-Dep.mat', sample), 'sample', 'kappa', 'dr', ...
      'settings', 'experiment', 'data', 'TInterval', 'NSteps', ...
      'L1', 'L2', 'Ld', 'Lh', 'N_e', 'dr_e', 'dr_p', 'k2D_e', 'k2D_p');
  else
    save(sprintf('water_T-Dep.mat'), 'sample', 'kappa', 'dr', ...
      'settings', 'experiment', 'data', 'TInterval', 'NSteps', ...
      'L1', 'L2', 'Ld', 'Lh', 'N_e', 'dr_e', 'dr_p', 'k2D_e', 'k2D_p');
  end
  clear p NTimeIntervals TInterval NSteps data dr_e dr_p 
  clear L1 L2 Ld Lh dr N_e k2D_e k2D_p kappa
end
clear sample settings experiment
       
