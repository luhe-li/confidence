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
    simD(k) = load(fileList(k).name);
    confCM = simD(k).conf_CM + confCM;
    locCM = simD(k).loc_CM + locCM;
end

%% merge files

loc_CM = zeros(3);
conf_CM = zeros(3);
for simd = 1:3
    for i_sample = 1:100
    % find the best fitting model for this sampled dataset
    [~, best_fitd] = min([saveLocModel{simd,1,i_sample}.minNLL, saveLocModel{simd,2,i_sample}.minNLL, saveLocModel{simd,3,i_sample}.minNLL]);
    loc_CM(simd, best_fitd) = loc_CM(simd, best_fitd) + 1;
    [~, best_fitd] = min([saveConfModel{simd,1,i_sample}.minNLL, saveConfModel{simd,2,i_sample}.minNLL, saveConfModel{simd,3,i_sample}.minNLL]);
    conf_CM(simd, best_fitd) = conf_CM(simd, best_fitd) + 1;
    end
end




figure(1); clf;
FM = round(100*CM/sum(CM(1,:)))/100;
t = imageTextMatrix(FM);
set(t(FM'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(CM);
set(t, 'fontsize', 22)
title(['count = ' num2str(count)]);
set(gca, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 28, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('simulated model')
