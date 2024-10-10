clear; close all; rng('shuffle');

%% knobs

n_sample = 1; % number of ground-truth samples to generate
n_run = 5;
reps = 3000; % number of trials per codnition
useCluster = false;
check_fake_data = false; % check the simulated data before fitting
i_model = 1;

%% model info

rng('Shuffle');
specifications = {'Unimodal, joint fit of localization and confidence'}; % official model names for plotting
folders = {'uni'}; % short names for folders
n_model = numel(specifications);
model_info = table((1:n_model)', specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});
curr_model_str = folders{i_model};

%% manage paths

[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
addpath(genpath(fullfile(project_dir, 'func')));
addpath(genpath(fullfile(git_dir, 'bads'))); % add optimization tool, here we use BADS for example
out_dir = fullfile(pwd, folders{i_model}); % output will be saved to the model folder
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

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
model.n_run = n_run; % number of initiation
model.modality = 2;

model.uni_sA = sA;
model.uni_sV = model.uni_sA;

model.aV = 1; % assume unbiased visual perceptual transformation for highly reliable stimulus
model.bV = 0;
model.sigma_motor = 1; % cm, should be estimated from pointing task data

%% set simulation parameterss
% 'aA','bA','\sigma_{A}','\sigma_{V]','\sigma_C','\sigma_{P}','\mu_P'
GT = [     1.1,    0.1,     10,    1,     5,      3,     0.5];
OPTIONS.TolMesh = 1e-4;

flnm = sprintf('%s-rep%i-%i', curr_model_str, min(reps), max(reps));

addpath(genpath(fullfile(pwd, curr_model_str)));
curr_func = str2func(['nll_' curr_model_str]);

for rr = 1:numel(reps)

    model.uni_nrep = reps(rr); % number of repetitions per condition
    t_sAV = sAV;

    for i_sample = 1:n_sample

        %% simulate fake data

        t_curr_func = curr_func;
        t_model = model;
        t_model.mode = 'predict';
        t_data = t_curr_func(GT, t_model);
        sim_data(rr, i_sample).data = t_data;
        sim_data(rr, i_sample).gt = GT;

        %% fit fake data

        t_model.mode = 'initiate';
        val = t_curr_func([], t_model, []);
        t_model.initVal = val;

        t_model.mode = 'optimize';
        llfun = @(x) t_curr_func(x, t_model, t_data);
        fprintf('[%s] Start parameter recover for model-%s, no. sample-%i \n', mfilename, curr_model_str, i_sample);

        test = llfun(GT);

        % fit the model multiple times with different initial values
        est_p = nan(t_model.n_run, val.num_param);
        nll = nan(1, t_model.n_run);
        parfor i  = 1:t_model.n_run
            [est_p(i,:), nll(i)] = bads(llfun,...
                val.init(i,:), val.lb, val.ub, val.plb, val.pub,[],OPTIONS);
        end

        % find the best fits across runs
        [min_nll, best_idx] = min(nll);
        best_p = est_p(best_idx, :);
        fits(rr, i_sample).est_p = est_p;
        fits(rr, i_sample).best_p = best_p;
        fits(rr, i_sample).nll = nll;
        fits(rr, i_sample).min_nll = min_nll;

    end

end

% save partial results
save(fullfile(out_dir, flnm), 'sim_data','fits');