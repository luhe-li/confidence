clear; close all; rng('shuffle');

%% model info
rng('Shuffle');
specifications = {'MA, full posterior, global','MA, full posterior, local','MA, gaussian posterior','MS','PM'};
folders = {'MA_optimal', 'MA_local', 'MA_gauss', 'MS','PM'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manege path
restoredefaultpath;
[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
% addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(pwd, 'bads')));
addpath(genpath(fullfile(git_dir, 'bads')));
out_dir = fullfile(pwd, mfilename);
if ~exist(outDir, 'dir'); mkdir(out_dir); end

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
delta_cm = unique(round(abs(sA - sV'),2));
delta_dva = rad2deg(atan(delta_cm / 2 / sitting_dist)) .* 2;

%% set simulation parameters

%         aA,     bA,  sigV1,   sigA,    sigP, sigConf,     pCC
GT = [     1,    0.1,     1,      10,      30,       1,     0.7];
n_para = length(GT);

% conditions
cue_label = {'Auditory post-cue','Visual post-cue'};
n_cue = numel(cue_label);
reps = [20, 30, 50, 70, 1000];

% simulate and fit for 100 times
n_sim = 100;

for fit_m = 1:numel(folders)

    for i_sample = 1:n_sim

        n_rep = reps(rr);
        curr_model_str = folders{fit_m};

        %% simulate fake data

        fixP.bi_sA = sAV(1,:);
        fixP.bi_sV = sAV(2,:);
        fixP.n_sA = numel(sA);
        fixP.sim_d = sim_d;
        fixP.bi_nrep = n_rep;
        fixP.screen_cm = screen_cm;
        fixP.center_axis = linspace(-screen_cm/2, screen_cm/2, 1e3);
        fixP.maxScore = 1;
        fixP.minScore = 0.01;
        fixP.elbow = screen_cm/4; % point goes to 0.01 when confidence range is 1/4 of screen
        fixP.dropRate = (fixP.maxScore - fixP.minScore)/fixP.elbow;

        % fixed parameters
        muP = 0;
        fixP.sigMotor = 3; % localization noise, in cm

        addpath(genpath(fullfile(project_dir, folders{fit_m})));
        curr_func = str2func(['sim_' folders{fit_m}]);

        [org_loc, org_conf] = curr_func(GT, fixP);
        fake_data(sim_m, i_sample).org_loc = org_loc;
        fake_data(sim_m, i_sample).org_conf = org_conf;

        %% fit fake data

        model.saveR = 0;
        model.mode = 'initialize';
        Val = curr_func([], model, []);
        model.initVal = Val;

        model.mode = 'optimize';

        printf('[%s] Start recover model-%s, recovery sample-%i \n', mfilename, folders{fit_m}, i_sample);


        % model prediction
    

    end
end