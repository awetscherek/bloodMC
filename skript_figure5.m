
% plots ADC0, ADC400 and Kapp data into four plots as in Fig. 5 in the
% manuscript. Uses the 's01_T-Dep_results.mat' files in the current
% directory. So if this script is executed in the 300k subdirectory, Fig. 5
% is reproduced.

clear
load('s01_T-Dep_results.mat')

figure(1)
ColOrd = get(gca,'ColorOrder');
hold off
plot(experiment.T_, Da_mp(:, 4), 'Color', ColOrd(1, :))
hold on
plot(experiment.T_, Dr_mp(:, 4), '+', 'Color', ColOrd(1, :))

figure(2)
hold off
plot(experiment.T_, Da_fc(:, 4), 'Color', ColOrd(1, :))
hold on
plot(experiment.T_, Dr_fc(:, 4), '+', 'Color', ColOrd(1, :))

figure(3)
hold off
plot(experiment.T_, Ka_mp(:, 4), 'x', 'Color', ColOrd(1, :))
hold on

figure(4)
hold off
plot(experiment.T_, Ka_fc(:, 4), 'x', 'Color', ColOrd(1, :))
hold on

sample = 1;
while (exist(sprintf('s%02d_T-Dep_results.mat', sample + 1), 'file'))  
  sample = sample + 1;
  load(sprintf('s%02d_T-Dep_results.mat', sample))
  
  figure(1)
  plot(experiment.T_, Da_mp(:, 4), 'Color', ColOrd(sample, :))
  plot(experiment.T_, Dr_mp(:, 4), '+', 'Color', ColOrd(sample, :))
  
  figure(2)
  plot(experiment.T_, Da_fc(:, 4), 'Color', ColOrd(sample, :))
  plot(experiment.T_, Dr_fc(:, 4), '+', 'Color', ColOrd(sample, :))
  
  figure(3)
  plot(experiment.T_, Ka_mp(:, 4), 'x', 'Color', ColOrd(sample, :))
  
  figure(4)
  plot(experiment.T_, Ka_fc(:, 4), 'x', 'Color', ColOrd(sample, :))   
end
