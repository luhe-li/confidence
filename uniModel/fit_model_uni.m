clear; close all; clc;

%% knobs
i_model = 1; % selected model to fit

%% specify models
rng('shuffle');

specifications = {'Unimodal, joint fit of localization and confidence'}; % official model names for plotting
folders = {'uni'}; % short names for folders
n_model = numel(specifications);
model_info = table((1:n_model)', specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
addpath(genpath(fullfile(project_dir, 'util')));
addpath(genpath(fullfile(git_dir, 'bads'))); % add optimization tool, here we use BADS for example
out_dir = fullfile(pwd, folders{i_model}); % output will be saved to the model folder
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% organize data

[data, exp_info] = org_resp('OY', {'A','V'}, 'uniLoc');

%% experiment info

% cm
speaker_cm = 65.5; % cm, left to center
sitting_dist = 113; % cm
screen_cm = 170; % cm

% pixel
screen_px = 1024; % pixel
px_axis = 1:screen_px; % add fence
pixel_per_cm = screen_px/screen_cm;

% dva
screen_dva = rad2deg(atan(screen_cm / 2 / sitting_dist)) .* 2;

% stimulus location in pixel, dva
aud_level = [5 7 10 12];
speaker_level = linspace(-speaker_cm, speaker_cm, 16);
sA = speaker_level(aud_level);
sA_dva = rad2deg(atan(sA / 2 / sitting_dist)) .* 2;
sV = sA; % assume audiovisual bias corrected
sAV = combvec(sA, sV);
model.delta_cm = unique(round(abs(sA - sV'),2));
model.delta_dva = rad2deg(atan(model.delta_cm / 2 / sitting_dist)) .* 2;

%% model setup

model.screen_cm = screen_cm;
model.center_axis = linspace(-screen_cm/2, screen_cm/2, 1e3);
model.maxScore = 1;
model.minScore = 0.01;
model.elbow = screen_cm/4; % point goes to 0.01 when confidence range is 1/4 of screen
model.dropRate = (model.maxScore - model.minScore)/model.elbow;
model.center_axis = linspace(-screen_cm/2, screen_cm/2, 1e3);
model.modality = 2;

model.uni_sA = sA;
model.uni_sV = model.uni_sA;
model.uni_nrep = exp_info.nRep; % number of repetitions per condition

model.n_run = 5; % number of fits for each model
model.aV = 1; % assume unbiased visual perceptual transformation for highly reliable stimulus
model.bV = 0;
model.sigma_motor = 1; % cm, should be estimated from pointing task data

%% fit model

fit_str = folders{i_model};
curr_model = str2func(['nll_' fit_str]);

model.mode = 'initiate';
val = curr_model([], model, []);
model.init_val = val;

model.mode = 'optimize';
llfun = @(x) curr_model(x, model, data);
fprintf('[%s] Start fitting model-%s\n', mfilename, fit_str);

% fit the model multiple times with different initial values
est_p = nan(model.n_run, val.num_param);
nll = nan(1, model.n_run);
for i  = 1:model.n_run
    [est_p(i,:), nll(i)] = bads(llfun,...
        val.init(i,:), val.lb, val.ub, val.plb, val.pub);
end

% find the best fits across runs
[min_nll, best_idx] = min(nll);
best_p = est_p(best_idx, :);
fits.best_p = best_p;
fits.min_nll = min_nll;

%% model prediction using the best-fitting parameters

model.mode = 'predict';
pred = curr_model(best_p, model, []);

%% save full results
fprintf('[%s] Model recovery done! Saving full results.\n', mfilename);
flnm = 'example_results';
save(fullfile(out_dir, flnm), 'data','model','fits','pred');
