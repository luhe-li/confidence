clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
results_dir           = fullfile(cur_dir, 'modelRecovery');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% load fitting results

num_rep = 30; % experimental repetition per condition
num_run = 10; % run of fits
num_sample = 100; % sample of GT
flnm = sprintf('RecoveryResults_*_rep%i_sample%i_run%i.mat', num_rep, num_sample, num_run);
fileList = dir(fullfile(results_dir,flnm));

confCM = zeros(3);
locCM = zeros(3);
for k = 1:length(fileList)
    simD(k) = load(fullfile(results_dir, fileList(k).name));
    confCM = simD(k).conf_CM + confCM;
    locCM = simD(k).loc_CM + locCM;
end

%% merge files

loc_CM = zeros(3);
conf_CM = zeros(3);
for simd = 1:3
    for i_sample = 1:100
    % find the best fitting model for this sampled dataset
    [~, best_fitd1] = min([simD(simd).saveLocModel{simd,1,i_sample}.minNLL, simD(simd).saveLocModel{simd,2,i_sample}.minNLL, simD(simd).saveLocModel{simd,3,i_sample}.minNLL]);
    loc_CM(simd, best_fitd1) = loc_CM(simd, best_fitd1) + 1;
    [~, best_fitd] = min([simD(simd).saveConfModel{simd,1,i_sample}.minNLL, simD(simd).saveConfModel{simd,2,i_sample}.minNLL, simD(simd).saveConfModel{simd,3,i_sample}.minNLL]);
    conf_CM(simd, best_fitd) = conf_CM(simd, best_fitd) + 1;
    end
end

%% plot

figure(1); clf;
t = imageTextMatrix(loc_CM);
set(t(loc_CM'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(loc_CM);
set(t, 'fontsize', 22)
title('Localization fit, sim = 100', 'fontsize',25)
set(gca, 'xtick', [1:3], 'ytick', [1:3], 'fontsize', 15, ...
    'xaxislocation', 'top', 'tickdir', 'out')
set(gca, 'xticklabel', {'heuristic', 'suboptimal', 'optimal'}, ...
    'yticklabel', {'heuristic', 'suboptimal', 'optimal'})
xlabel('Fit model', 'fontsize', 20)
ylabel('Simulated model', 'fontsize', 20)
saveas(gca, fullfile(out_dir, sprintf('locCM_rep-%i_sample-%i', num_rep, num_sample)), 'png')

figure(2); clf;
t = imageTextMatrix(conf_CM);
set(t(conf_CM'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(conf_CM);
set(t, 'fontsize', 22)
title('Localization & confidence fit, sim = 100', 'fontsize',25)
set(gca, 'xtick', [1:3], 'ytick', [1:3], 'fontsize', 15, ...
    'xaxislocation', 'top', 'tickdir', 'out')
set(gca, 'xticklabel', {'heuristic', 'suboptimal', 'optimal'}, ...
    'yticklabel', {'heuristic', 'suboptimal', 'optimal'})
xlabel('Fit model', 'fontsize', 20)
ylabel('Simulated model', 'fontsize', 20)
saveas(gca, fullfile(out_dir, sprintf('confCM_rep-%i_sample-%i', num_rep, num_sample)), 'png')

