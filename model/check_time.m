clear; close all; rng('shuffle');

%% knobs

n_sample = 1; % number of ground-truth samples to generate
n_run = 1;
n_rep = 30; % number of trials per condition
useCluster = false;
check_fake_data = false; % check the simulated data before fitting

%% model info

rng('Shuffle');
specifications = {'MA, full posterior, global','MA, full posterior, local','MA, gaussian posterior','MS','PM'};
folders = {'MA_optimal', 'MA_local', 'MA_gauss', 'MS','PM'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% set environment

if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster= false;
end

% job = model
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);
        model_slc  = hpc_job_number;

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false
        % for local debug
        numCores = feature('numcores');
        model_slc = 5;
end

%% manage path

restoredefaultpath;
[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
% addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(pwd, 'util')));
addpath(genpath(fullfile(git_dir, 'bads')));
out_dir = fullfile(pwd, mfilename);
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
model.n_run = n_run;%3; % number of initiation
model.num_SD = 5;
model.numBins_A = 15;
model.numBins_V = 15;
model.modality = 2;

% fixed params
model.sigMotor = 3; % localization noise, in cm
model.muP = 0;

% conditions
model.cue_label = {'Auditory post-cue','Visual post-cue'};
model.n_cue = numel(model.cue_label);

%% set simulation parameters

%         aA,     bA,  sigV1,   sigA,    sigP, sigConf,     pCC
GT = [     1,    0.1,     1,      10,      10,       1,     0.7];
OPTIONS.TolMesh = 1e-4;

time_results = zeros(numel(folders), 1);
sim_data = struct();
fits = struct();

parfor mm = 1:numel(folders)
    t_sAV = sAV;
    t_model = model;
    t_model.bi_sA = t_sAV(1,:); % 1x16
    t_model.bi_sV = t_sAV(2,:);
    t_model.sA = sA; %1x4
    t_model.sV = sA;
    t_model.n_sA = numel(sA);
    t_model.bi_nrep = n_rep;

    addpath(genpath(fullfile(pwd, folders{mm})));
    curr_func = str2func(['NLL_' folders{mm}]);

    for rr = 1
        for i_sample = 1:n_sample
            tic; % Start timing
            
            %% simulate fake data

            t_curr_func = curr_func;
            t_model.mode = 'predict';
            t_data = t_curr_func(GT, t_model);
            sim_data(mm, rr, i_sample).data = t_data;
            sim_data(mm, rr, i_sample).gt = GT;

            %% fit fake data

            t_model.saveR = 0;
            t_model.mode = 'initiate';
            val = t_curr_func([], t_model, []);
            t_model.initVal = val;

            t_model.mode = 'optimize';
            llfun = @(x) t_curr_func(x, t_model, t_data);
            fprintf('[%s] Start parameter recover for model-%s, no. sample-%i \n', mfilename, folders{mm}, i_sample);

            test = llfun(GT);

            % fit the model multiple times with different initial values
            est_p = nan(model.n_run, val.num_param);
            nll = nan(1, val.num_param);
            for i  = 1:model.n_run
                t_val = val;
                [est_p(i,:), nll(i)] = bads(llfun,...
                    t_val.init(i,:), t_val.lb, t_val.ub, t_val.plb, t_val.pub,[],OPTIONS);
            end

            % find the best fits across runs
            [min_nll, best_idx] = min(nll);
            best_p = est_p(best_idx, :);
            fits(mm, rr, i_sample).est_p = est_p;
            fits(mm, rr, i_sample).best_p = best_p;
            fits(mm, rr, i_sample).nll = nll;
            fits(mm, rr, i_sample).min_nll = min_nll;

            time_results(mm) = toc; % End timing and store result
        end
    end
end

%% save full results

fprintf('[%s] Timing check done! Saving full results.\n', mfilename);
flnm = sprintf('%s-rep%i-%i', curr_model_str, n_rep, n_rep);
save(fullfile(out_dir, flnm), 'sim_data','fits', 'time_results');

% Display timing results
for mm = 1:numel(folders)
    fprintf('Model %s took %.2f seconds to fit.\n', folders{mm}, time_results(mm));
end