clear; close all; rng('shuffle');

sub_slc = 12;
ses_slc = 1; % bimodal sessions

models = {'Heuristic','Suboptimal','Optimal'};
d = 3;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
result_dir            = fullfile(cur_dir, 'modelFit');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(result_dir))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;
%     251, 154, 153]./255; % light red

%% predict data

% load best fitting parameters
flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
files = dir(fullfile(result_dir, flnm));
load(files(end).name);

% use the best-fitting parameter and a selected model
p = saveConfModel{d}.bestP;
model.mode = 'predict';
model.model_slc             = d;
model.strategy_conf         = models{d};
pred = nllLocConf(p, model, data);

%% analyze data prediction 

pred.biExpInfo = data.biExpInfo;

[pred_mean_conf, pred_std_mean_conf, pred_uni_pconf, ~,...
    pred_mean_ve, pred_std_ve, ~] = analyze_data(sub_slc, ses_slc, pred);

%% analyze real data

[mean_conf, std_mean_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc);


