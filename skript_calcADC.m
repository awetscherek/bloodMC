       
%%
clear

% script assumes that results from T-Dep simulations are available:
% due to memory restrictions each sample file probably needs to be processed
% successively.
sample =    1;
load(sprintf('s%02d_T-Dep.mat', sample), 'data', 'TInterval', 'NSteps', ...
    'settings', 'experiment', 'L1', 'L2', 'Ld', 'Lh', 'N_e', 'kappa');

%sample =    01;
%load(sprintf('water_T-Dep.mat'), 'data', 'TInterval', 'NSteps', ...
%    'settings', 'experiment', 'L1', 'L2', 'Ld', 'Lh', 'N_e', 'kappa');

% extract everything except positions at the end of intervals:
time_factor = repmat(TInterval.' ./ NSteps.', size(data, 1), 1);
x0     = data(:, 1);                                   % starting pos
y0     = data(:, 2);                           
z0     = data(:, 3);
xcm     = data(:, 4:9:size(data, 2));                  % mean pos in each interval
ycm     = data(:, 5:9:size(data, 2));
zcm     = data(:, 6:9:size(data, 2));
t_in    = data(:, 7:9:size(data, 2)) .* time_factor;   % time spent withing RBC
t_out   = data(:, 8:9:size(data, 2)) .* time_factor;   % time spent in Plasma
t_first = data(:, 9:9:size(data, 2)) .* time_factor;   % first change of compartments
clear data time_factor NSteps

fprintf('\n--------------------------------------------------------\n');

fprintf(['\nSample No.%2d, geometric HCT = % 3.3f %%, ' ...
      'MCV = % 3.1f fL\n'], sample, 100 * (Ld.^2 .* Lh * pi / 4) ...
      ./ (L1.^2 .* L2), pi * Ld.^2 .* Lh / 4e-18); % nominal HCT
  
fprintf('\nHCT from residence times:\n');
for j = 1:numel(TInterval)   
  fprintf('Interval %2d: % 3.3f %%\n', j, settings.c_p * mean(t_in(:, j)) ...
    / TInterval(j) / (settings.c_e + (settings.c_p - settings.c_e) ...
    * mean(t_in(:, j)) / TInterval(j)) * 100);     % HCT based on trajectories / concentrations
end
clear j
fprintf('Average HCT from residence times = % 3.3f %%\n', settings.c_p ...
    * sum(mean(t_in, 1)) * 100 / sum(TInterval) / (settings.c_e ...
    + (settings.c_p-settings.c_e) * sum(mean(t_in, 1)) / sum(TInterval)));

outside = ((x0-L1/2).^2 + (y0-L1/2).^2 > Ld^2/4) | ((z0-L2/2).^2 > Lh^2/4);
fprintf('\nFound %d particles starting in Plasma ... ', sum(outside));
if (size(x0, 1) - sum(outside) == N_e)
    fprintf('as expected\n');
else
    fprintf('\n');
    warning('expected %d particles', size(x0, 1) - N_e);
end
clear outside x0 y0 z0 L1 L2

first_out = zeros(N_e, 1);
for p = 1:N_e   
    s = 1;
    while ((s <= size(t_first, 2)) && (t_first(p, s) == 0))
        s = s + 1;     
    end
    if (s <= size(t_first, 2))
        first_out(p, 1) = sum(t_in(p, 1:(s-1))) + t_first(p, s);
    else
        first_out(p, 1) = -1;
    end
end

first_in = zeros(size(t_first, 1) - N_e, 1);
for p = (N_e+1):size(t_first, 1)
    s = 1;
    while ((s <= size(t_first, 2)) && (t_first(p, s) == 0))
        s = s + 1;     
    end
    if (s <= size(t_first, 2))
        first_in(p - N_e, 1) = sum(t_out(p, 1:(s-1))) - t_first(p, s);
    else
        first_in(p - N_e, 1) = -1;
    end  
end
clear p s t_first

arrival_times    = [0; sort(first_in(first_in  > 0))].';
remaining_plasma = numel(first_in ) - (0:numel(find(first_in  > 0)));

departure_times  = [0; sort(first_out(first_out > 0))].';
remaining_RBC    = numel(first_out) - (0:numel(find(first_out > 0)));
clear first_in first_out

R_p = fminsearch(@(x) norm(remaining_plasma - ...
    remaining_plasma(1) * exp(-arrival_times   * x)), 0); 

R_e = fminsearch(@(x) norm(remaining_RBC    - ...
    remaining_RBC(1)    * exp(-departure_times * x)), 0); 

plot(arrival_times,   remaining_plasma, '.');
hold on
plot(departure_times, remaining_RBC   , '.');

set(gca, 'ColorOrderIndex', 1)

pnts = 1000;
plot((0:pnts) * max(arrival_times) / pnts, ...
    remaining_plasma(1) * exp(-(0:pnts) * max(arrival_times) * R_p / pnts));
plot((0:pnts) * max(departure_times) / pnts, ...
    remaining_RBC(1) * exp(-(0:pnts) * max(departure_times) * R_e / pnts));
hold off
legend('Particles in Plasma', 'Particles in RBC', ...
  sprintf('fit R_p = % 3.3f Hz',R_p), sprintf('fit R_e = % 3.3f Hz',R_e));
clear arrival_times departure_times pnts remaining_RBC remaining_plasma

fprintf('Expected  intra-cellular pre-exchange lifetime = %f ms\n', ...
  Ld * Lh / (2 * Ld + 4 * Lh) / kappa * 1000);
clear Ld Lh kappa

fprintf('Simulated intra-cellular pre-exchange lifetime = %f ms\n', ... 
  1000 / R_e);
clear R_e

fprintf('Simulated extra-cellular pre-exchange lifetime = %f ms\n', ... 
  1000 / R_p);
clear R_p

% variables to store for each interval the displacement as a sum of the mean 
% positions weighted by the gradient sign. This property is proportional to
% the acquired phase.
stejskal_delta_x = zeros(size(xcm, 1), numel(experiment.T_));
stejskal_delta_y = zeros(size(ycm, 1), numel(experiment.T_));
stejskal_delta_z = zeros(size(zcm, 1), numel(experiment.T_));

fc_plus_delta_x = zeros(size(xcm, 1), numel(experiment.T_));
fc_plus_delta_y = zeros(size(ycm, 1), numel(experiment.T_));
fc_plus_delta_z = zeros(size(zcm, 1), numel(experiment.T_));

tweight = ones(size(xcm, 1), 1) * TInterval.'; % weight each interval by its duration

for T = experiment.T_
    
   % we only need to look at the intervals within the diffusion time, so we calculate those...     
   [~, TInterval2] = ...
      calcIntervalDurations(experiment.TE * 1e-3, T * 1e-3);
  
   % ... and determine where the first gradient lobe starts:
   fc1_start = 1;
   while (TInterval2(1) - sum(TInterval(1:fc1_start)) > 1e-10)
       fc1_start = fc1_start + 1;
   end
   fc1_start = fc1_start + 1;
      
   % ... and where it ends:
   fc1_end = fc1_start;
   while (TInterval2(2) - sum(TInterval(fc1_start:fc1_end)) > 1e-10)
       fc1_end = fc1_end + 1;
   end
    
   % ... and the second one:
   fc2_start = fc1_end + 1;
   fc2_end   = fc2_start;
   while (TInterval2(3) - sum(TInterval(fc2_start:fc2_end)) > 1e-10)
       fc2_end = fc2_end + 1;
   end
   
   % third one:
   fc3_start = fc2_end + 1;
   fc3_end   = fc3_start;
   while (TInterval2(4) - sum(TInterval(fc3_start:fc3_end)) > 1e-10)
       fc3_end = fc3_end + 1;
   end
   
   % last one:
   fc4_start = fc3_end + 1;
   fc4_end   = fc4_start;
   while (TInterval2(5) - sum(TInterval(fc4_start:fc4_end)) > 1e-10)
       fc4_end = fc4_end + 1;
   end
   
   % noe we can calculate the weighted sum. For monopolar the first two
   % intervals of the flow-compensated can be merged:
   stejskal_delta_x(:, T == experiment.T_) = 1000 * (...
     sum(xcm(:,fc1_start:fc2_end) .* tweight(:,fc1_start:fc2_end), 2) - ...
     sum(xcm(:,fc3_start:fc4_end) .* tweight(:,fc3_start:fc4_end), 2)) / T;
 
   stejskal_delta_y(:, T == experiment.T_) = 1000 * (...
     sum(ycm(:,fc1_start:fc2_end) .* tweight(:,fc1_start:fc2_end), 2) - ...
     sum(ycm(:,fc3_start:fc4_end) .* tweight(:,fc3_start:fc4_end), 2)) / T;
 
   stejskal_delta_z(:, T == experiment.T_) = 1000 * (...
     sum(zcm(:,fc1_start:fc2_end) .* tweight(:,fc1_start:fc2_end), 2) - ...
     sum(zcm(:,fc3_start:fc4_end) .* tweight(:,fc3_start:fc4_end), 2)) / T;
 
   % for flow-compensated the 2nd and 4th lobe have negative polarity:
   fc_plus_delta_x (:, T == experiment.T_) = 1000 * (...
     sum(xcm(:,fc1_start:fc1_end) .* tweight(:,fc1_start:fc1_end), 2) - ...
     sum(xcm(:,fc2_start:fc2_end) .* tweight(:,fc2_start:fc2_end), 2) + ...
     sum(xcm(:,fc3_start:fc3_end) .* tweight(:,fc3_start:fc3_end), 2) - ...
     sum(xcm(:,fc4_start:fc4_end) .* tweight(:,fc4_start:fc4_end), 2)) / T;
 
   fc_plus_delta_y (:, T == experiment.T_) = 1000 * (...
     sum(ycm(:,fc1_start:fc1_end) .* tweight(:,fc1_start:fc1_end), 2) - ...
     sum(ycm(:,fc2_start:fc2_end) .* tweight(:,fc2_start:fc2_end), 2) + ...
     sum(ycm(:,fc3_start:fc3_end) .* tweight(:,fc3_start:fc3_end), 2) - ...
     sum(ycm(:,fc4_start:fc4_end) .* tweight(:,fc4_start:fc4_end), 2)) / T;
 
   fc_plus_delta_z (:, T == experiment.T_) = 1000 * (...
     sum(zcm(:,fc1_start:fc1_end) .* tweight(:,fc1_start:fc1_end), 2) - ...
     sum(zcm(:,fc2_start:fc2_end) .* tweight(:,fc2_start:fc2_end), 2) + ...
     sum(zcm(:,fc3_start:fc3_end) .* tweight(:,fc3_start:fc3_end), 2) - ...
     sum(zcm(:,fc4_start:fc4_end) .* tweight(:,fc4_start:fc4_end), 2)) / T;
    
end
clear fc1_start fc1_end fc2_start fc2_end fc3_start fc3_end 
clear fc4_start fc4_end xcm ycm zcm tweight TInterval TInterval2

t_in  = sum(t_in,  2);
t_out = sum(t_out, 2);

% initialize variables to store the results:
Da_mp    = zeros(numel(experiment.T_), 4);
Da_fc    = zeros(numel(experiment.T_), 4);
Dr_mp    = zeros(numel(experiment.T_), 4);
Dr_fc    = zeros(numel(experiment.T_), 4);
Ka_mp    = zeros(numel(experiment.T_), 4);
Ka_fc    = zeros(numel(experiment.T_), 4);
mu_mp    = zeros(numel(experiment.T_), 4);
mu_fc    = zeros(numel(experiment.T_), 4);
sigma2_mp = zeros(numel(experiment.T_), 4);
sigma2_fc = zeros(numel(experiment.T_), 4);

% calculate centers of masses along many directions to get trace-weighting
NDirections =  10000; 
kk_theta = rand(1, NDirections) * 2 - 1;
kk_phi   = rand(1, NDirections) * 2 * pi;

% from appendix of MRM paper (magic_factor = b/gamma^2/g^2/T^3)
magic_factor_st = 1/12;
magic_factor_fc = sqrt(2)/8 - 1/6;

b = 400 * 1e6; % SI units
gradampl_larmor_T_st = sqrt(b / magic_factor_st ./ experiment.T_ / 1e-3);
gradampl_larmor_T_fc = sqrt(b / magic_factor_fc ./ experiment.T_ / 1e-3);


% assuming that there is no T2 relaxation:
R2Plasma = 0;
R2RBC    = 0;
w = repmat(exp(-t_in * R2RBC -t_out * R2Plasma), 1, NDirections);

%
% due to memory constraints first stejskal:

for t_index = 1:size(stejskal_delta_x, 2)
    
  % trace-weighted:
  stejskal_delta_t = zeros(size(stejskal_delta_x, 1), NDirections);
  
  for split_index = 1:100 
      
    split_range = (1:(size(stejskal_delta_x, 1) / 100)) + ...
        (size(stejskal_delta_x, 1) / 100) * (split_index - 1);
    
    stejskal_delta_t(split_range, :) = ...
      stejskal_delta_x(split_range, t_index) * ...
        (cos(kk_phi).* sqrt(1-kk_theta.^2)) + ...
      stejskal_delta_y(split_range, t_index) * ...
        (sin(kk_phi).* sqrt(1-kk_theta.^2)) + ...
      stejskal_delta_z(split_range, t_index) * kk_theta;
  
  end

   T = experiment.T_(t_index);
   sumw1 = sum(w(:, 1));
   sumw  = sum(w(:));
   
   Dr_mp(t_index, 1) = -log(abs(sum(exp(1i * stejskal_delta_x(:, t_index) ...
    * gradampl_larmor_T_st(t_index)).* w(:, 1)) / sumw1)) / b;
   Dr_mp(t_index, 2) = -log(abs(sum(exp(1i * stejskal_delta_y(:, t_index) ...
    * gradampl_larmor_T_st(t_index)).* w(:, 1)) / sumw1)) / b;
   Dr_mp(t_index, 3) = -log(abs(sum(exp(1i * stejskal_delta_z(:, t_index) ...
    * gradampl_larmor_T_st(t_index)).* w(:, 1)) / sumw1)) / b;
  
   Dr_mp(t_index, 4) = -log(abs(sum(exp(1i * stejskal_delta_t(:) ...
    * gradampl_larmor_T_st(t_index)) .* w(:)) / sumw)) / b;
 
   mu_mp(t_index, 1) = sum(w(:, 1) .* stejskal_delta_x(:, t_index)) / sumw1;
   mu_mp(t_index, 2) = sum(w(:, 1) .* stejskal_delta_y(:, t_index)) / sumw1;
   mu_mp(t_index, 3) = sum(w(:, 1) .* stejskal_delta_z(:, t_index)) / sumw1;
   mu_mp(t_index, 4) = sum(w(:)    .* stejskal_delta_t(:)) / sumw;
   
   sigma2_mp(t_index, 1) = sum(w(:, 1) .* (stejskal_delta_x(:, t_index) - mu_mp(t_index, 1)).^2) / sumw1;
   sigma2_mp(t_index, 2) = sum(w(:, 1) .* (stejskal_delta_y(:, t_index) - mu_mp(t_index, 2)).^2) / sumw1; 
   sigma2_mp(t_index, 3) = sum(w(:, 1) .* (stejskal_delta_z(:, t_index) - mu_mp(t_index, 3)).^2) / sumw1;
   sigma2_mp(t_index, 4) = sum(w(:)    .* (stejskal_delta_t(:)          - mu_mp(t_index, 4)).^2) / sumw;
   
   Da_mp(t_index, 1) = 1000 / T / magic_factor_st * sigma2_mp(t_index, 1) / 2;
   Da_mp(t_index, 2) = 1000 / T / magic_factor_st * sigma2_mp(t_index, 2) / 2;
   Da_mp(t_index, 3) = 1000 / T / magic_factor_st * sigma2_mp(t_index, 3) / 2;
   Da_mp(t_index, 4) = 1000 / T / magic_factor_st * sigma2_mp(t_index, 4) / 2;
   
   Ka_mp(t_index, 1) = sum(w(:, 1) .* (stejskal_delta_x(:, t_index) - mu_mp(t_index, 1)).^4) / sumw1 / sigma2_mp(t_index, 1).^2 - 3;
   Ka_mp(t_index, 2) = sum(w(:, 1) .* (stejskal_delta_y(:, t_index) - mu_mp(t_index, 2)).^4) / sumw1 / sigma2_mp(t_index, 2).^2 - 3;
   Ka_mp(t_index, 3) = sum(w(:, 1) .* (stejskal_delta_z(:, t_index) - mu_mp(t_index, 3)).^4) / sumw1 / sigma2_mp(t_index, 3).^2 - 3;
   Ka_mp(t_index, 4) = sum(w(:)    .* (stejskal_delta_t(:)          - mu_mp(t_index, 4)).^4) / sumw  / sigma2_mp(t_index, 4).^2 - 3;
 
end

clear stejskal_delta_t

% then the same for fc
for t_index = 1:size(fc_plus_delta_x, 2)
    
  fc_plus_delta_t = zeros(size(fc_plus_delta_x, 1), NDirections);
  
  for split_index = 1:100 
      
    split_range = (1:(size(fc_plus_delta_x, 1) / 100)) + ...
        (size(fc_plus_delta_x, 1) / 100) * (split_index - 1);
    
    fc_plus_delta_t(split_range, :) = ...
      fc_plus_delta_x(split_range, t_index) * ...
        (cos(kk_phi).* sqrt(1-kk_theta.^2)) + ...
      fc_plus_delta_y(split_range, t_index) * ...
        (sin(kk_phi).* sqrt(1-kk_theta.^2)) + ...
      fc_plus_delta_z(split_range, t_index) * kk_theta;
  
  end

   T = experiment.T_(t_index);
   sumw1 = sum(w(:, 1));
   sumw  = sum(w(:));
   
   Dr_fc(t_index, 1) = -log(abs(sum(exp(1i * fc_plus_delta_x(:, t_index) ...
    * gradampl_larmor_T_fc(t_index)).* w(:, 1)) / sumw1)) / b;
   Dr_fc(t_index, 2) = -log(abs(sum(exp(1i * fc_plus_delta_y(:, t_index) ...
    * gradampl_larmor_T_fc(t_index)).* w(:, 1)) / sumw1)) / b;
   Dr_fc(t_index, 3) = -log(abs(sum(exp(1i * fc_plus_delta_z(:, t_index) ...
    * gradampl_larmor_T_fc(t_index)).* w(:, 1)) / sumw1)) / b;
  
   Dr_fc(t_index, 4) = -log(abs(sum(exp(1i * fc_plus_delta_t(:) ...
    * gradampl_larmor_T_fc(t_index)) .* w(:)) / sumw)) / b;
 
   mu_fc(t_index, 1) = sum(w(:, 1) .* fc_plus_delta_x(:, t_index)) / sumw1;
   mu_fc(t_index, 2) = sum(w(:, 1) .* fc_plus_delta_y(:, t_index)) / sumw1;
   mu_fc(t_index, 3) = sum(w(:, 1) .* fc_plus_delta_z(:, t_index)) / sumw1;
   mu_fc(t_index, 4) = sum(w(:)    .* fc_plus_delta_t(:)         ) / sumw;
   
   sigma2_fc(t_index, 1) = sum(w(:, 1) .* (fc_plus_delta_x(:, t_index) - mu_fc(t_index, 1)).^2) / sumw1;
   sigma2_fc(t_index, 2) = sum(w(:, 1) .* (fc_plus_delta_y(:, t_index) - mu_fc(t_index, 2)).^2) / sumw1; 
   sigma2_fc(t_index, 3) = sum(w(:, 1) .* (fc_plus_delta_z(:, t_index) - mu_fc(t_index, 3)).^2) / sumw1;
   sigma2_fc(t_index, 4) = sum(w(:)    .* (fc_plus_delta_t(:)          - mu_fc(t_index, 4)).^2) / sumw;
   
   Da_fc(t_index, 1) = 1000 / T / magic_factor_fc * sigma2_fc(t_index, 1) / 2;
   Da_fc(t_index, 2) = 1000 / T / magic_factor_fc * sigma2_fc(t_index, 2) / 2;
   Da_fc(t_index, 3) = 1000 / T / magic_factor_fc * sigma2_fc(t_index, 3) / 2;
   Da_fc(t_index, 4) = 1000 / T / magic_factor_fc * sigma2_fc(t_index, 4) / 2;
   
   Ka_fc(t_index, 1) = sum(w(:, 1) .* (fc_plus_delta_x(:, t_index) - mu_fc(t_index, 1)).^4) / sumw1 / sigma2_fc(t_index, 1).^2 - 3;
   Ka_fc(t_index, 2) = sum(w(:, 1) .* (fc_plus_delta_y(:, t_index) - mu_fc(t_index, 2)).^4) / sumw1 / sigma2_fc(t_index, 2).^2 - 3;
   Ka_fc(t_index, 3) = sum(w(:, 1) .* (fc_plus_delta_z(:, t_index) - mu_fc(t_index, 3)).^4) / sumw1 / sigma2_fc(t_index, 3).^2 - 3;
   Ka_fc(t_index, 4) = sum(w(:)    .* (fc_plus_delta_t(:)          - mu_fc(t_index, 4)).^4) / sumw  / sigma2_fc(t_index, 4).^2 - 3;
 
end

sigma_mp = sqrt(sigma2_mp);
sigma_fc = sqrt(sigma2_fc);

clear kk_theta kk_phi t_index sigma2_fc sigma2_mp


save(sprintf('s%02d_T-Dep_results.mat', sample), 'Ka_fc', 'Ka_mp', 'Da_fc', ...
    'Da_mp', 'Dr_fc', 'Dr_mp', 'mu_fc', 'mu_mp', 'sigma_mp', ...
    'sigma_fc', 'experiment', 'settings', 'b');

%save(sprintf('water_T-Dep_results.mat'), 'Ka_fc', 'Ka_mp', 'Da_fc', ...
%    'Da_mp', 'Dr_fc', 'Dr_mp', 'mu_fc', 'mu_mp', 'sigma2_mp', ...
%    'sigma2_fc', 'experiment', 'settings', 'b');

%% in case one would like to fuse two phase data files into one:

clear

input1 = '50k/s%02d_T-Dep.mat';
input2 = '250k/s%02d_T-Dep.mat';
output = '300k/s%02d_T-Dep.mat';

for sample = 1:10
   
    data1 = load(sprintf(input1, sample));
    data2 = load(sprintf(input2, sample));
    
    data3 = data1;  
    data3.N_e = data1.N_e + data2.N_e;
    data3.settings.NRep = data1.settings.NRep + data2.settings.NRep;
    
    data3.data = zeros(data3.settings.NRep, size(data3.data, 2));
    
    data3.data(1:data1.N_e, :)               = data1.data(1:data1.N_e, :);
    data3.data(data1.N_e + (1:data2.N_e), :) = data2.data(1:data2.N_e, :);
    data3.data(data1.N_e + data2.N_e + ...
        (1:(data1.settings.NRep - data1.N_e)), :) ...
        = data1.data(data1.N_e + (1:(data1.settings.NRep - data1.N_e)), :);
    data3.data(data1.settings.NRep + data2.N_e + ...
        (1:(data2.settings.NRep - data2.N_e)), :) ...
        = data2.data(data2.N_e + (1:(data2.settings.NRep - data2.N_e)), :);
    
    [pathstr, ~, ~] = fileparts(sprintf(output, sample));
    if (~exist(pathstr, 'dir'))
        mkdir(pathstr);
    end
    
    save(sprintf(output, sample), '-struct', 'data3');
end

clear data1 data2 data3 input1 input2 output sample