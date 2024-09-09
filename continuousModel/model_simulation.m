
clear; close all; rng('shuffle');

%% set environment

sim_d = 1;
useCluster = false;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

% job/sample = 100, core/run = 3, fit once
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);


        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores);
        end

    case false
        numCores = feature('numcores');

end

%% manage paths

restoredefaultpath;
[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'func')));
addpath(genpath(fullfile(git_dir, 'vbmc')));
outDir = fullfile(pwd, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end

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
aud_level = 5:12;%[5 8 9 12];
speaker_level = linspace(-speaker_cm, speaker_cm, 16);
sA = speaker_level(aud_level);
sA_dva = rad2deg(atan(sA / 2 / sitting_dist)) .* 2;
sV = sA; % assume audiovisual bias corrected
sAV = combvec(sA, sV);
delta_cm = unique(round(abs(sA - sV'),2));
delta_dva = rad2deg(atan(delta_cm / 2 / sitting_dist)) .* 2;

%% set simulation parameters

%         aA,     bA,  sigV1,   sigA,   sigV2,  sigP, sigConf,   pCC
GT = [     1,    0.1,     1,     10,     7,     30,      1,    0.7];
n_para = length(GT);

% main models
model_names = {'Optimal','Model average','Model selection','Heuristic'};
folders = {'optimal','MA','MS','heuristic'};
n_model = numel(model_names);

% conditions
cue_label = {'Auditory post-cue','Visual post-cue'};
n_cue = numel(cue_label);
rel_label = {'High visual reliability','Low visual reliability'};
n_rep = 30; % repitition per condition

%% simulate fake data

fixP.bi_sA = sAV(1,:);
fixP.bi_sV = sAV(2,:);
fixP.n_sA = numel(sA);
fixP.sim_d = sim_d;
fixP.bi_nrep = n_rep;
fixP.screen_cm = screen_cm;
% fixP.left_axis = linspace(1, screen_cm, 1e3);
fixP.center_axis = linspace(-screen_cm/2, screen_cm/2, 1e3);
fixP.maxScore = 1;
fixP.minScore = 0.01;
fixP.elbow = 170/3;
fixP.dropRate = (fixP.maxScore - fixP.minScore)/fixP.elbow;

aA                    = GT(1);
bA                    = GT(2);
sigV1                 = GT(3);
sigA                  = GT(4);
sigV2                 = GT(5);
sigP                  = GT(6);
sigConf               = GT(7);
pCommon               = GT(8);

muP = 0;
fixP.sigMotor = 1.36; % in deg, measured from first four participants

% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% n_rep
% [org_loc, org_conf] = deal(NaN(n_sa,n_sa, n_cue, n_rel, n_rep));
curr_func = str2func(['sim_' folders{sim_d}]);

[org_loc, org_conf] = curr_func(...
    aA, bA, sigA, sigV1, sigV2, muP, sigP, sigConf, pCommon, fixP);

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% check localization data

uni_loc = zeros(size(org_loc));

loc_a = repmat((sA*aA+bA)',[1,numel(sV)]);
loc_v = repmat(sV,[numel(sA),1]);

uni_loc(:,:,1,1:2,:) = repmat(loc_a, [1, 1, 1, 2, n_rep]);
uni_loc(:,:,2,1:2,:) = repmat(loc_v, [1, 1, 1, 2, n_rep]);

% loc at uni minus loc at bi
ve =  mean(org_loc,5) - mean(uni_loc, 5);

% diff x cue x reliability
% assume participants localized perfectly in the unisensory
% condition
[ve_by_raw_diff, raw_diff] = org_by_raw_diffs_4D(ve, sA);

% organize confidence data: {diff} cue x reliability x rep
[conf_by_diff, diff] = org_by_diffs(org_conf, sA);

save(sprintf('sim_%s.mat',folders{sim_d}),'org_loc','org_conf','ve_by_raw_diff','raw_diff','conf_by_diff','diff','model_names','sim_d','cue_label','rel_label','n_rep','fixP');

%% double check here

figure; hold on
t = tiledlayout(1, 2);
title(t,sprintf('%s, rep: %i', model_names{sim_d}, n_rep))
xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
ylabel(t, 'Shift of localization');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:n_cue
    nexttile
    title(cue_label{cue})
    axis equal
    hold on

    for rel = 1: 2

        i_ve = squeeze(ve_by_raw_diff(:, cue, rel));
        plot(raw_diff, i_ve, 'Color',clt(rel+1,:))

    end
    xlim([min(raw_diff), max(raw_diff)])
    xticks(raw_diff)
    yline(0,'--')
    if cue == 1
        plot(raw_diff, raw_diff,'k--')
    else
        plot(raw_diff, -raw_diff,'k--')
    end
end
