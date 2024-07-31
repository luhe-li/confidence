
clear; close all; rng('shuffle');

%% set environment

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
speaker_cm          = 65.5; % cm, left to center
sitting_dist        = 113; % cm
screen_cm           = 170; % cm

% pixel
screen_px           = 1024; % pixel
pixel_per_cm = screen_px/screen_cm;

% dva
screen_dva          = rad2deg(atan(screen_cm / 2 / sitting_dist)) .* 2;

% stimulus location in pixel
aud_level           = [5 8 9 12];
speaker_level = linspace(-speaker_cm, speaker_cm, 16);
sA = round(speaker_level(aud_level) * pixel_per_cm);
sV = sA; % assume audiovisual bias corrected
sAV = combvec(sA, sV);
n_s = size(sAV,2);% number of stimulus conditions


%% 

sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA; % assume corrected by audiovisual bias
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of stimulus conditions

diff                 = zeros(length(sA), length(sV));
for ii                 = 1:length(sA)
    for j                 = 1:length(sV)
        diff(ii, j)           = sA(ii) - sV(j);
    end
end
all_diffs             = unique(abs(diff))';

%% set simulation parameters

%         aA,     bA, sigV1,  sigA,   sigV2,  sigP,  sigConf, sig pCC
GT = [   1.1,   0.1,    1,    3,    4,    8,   0.6,   4,    0.5];
n_para    = length(GT);

ds_loc                = {'Model averaging'};
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
n_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
n_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

sim_d = 3;
n_rep = 30;

%% simulate fake data

aA                    = GT(1);
bA                    = GT(2);
sigV1                 = GT(3);
sigA                  = GT(4);
sigV2                 = GT(5);
sigVs                 = [sigV1, sigV2];
sigP                  = GT(6);
sigConf               = GT(7);
pCommon               = GT(8);

lapse                 = 0.02;
muP                   = 0;
sigLoc                = 4.4;

% num_s: stimulus location combination
% 3 confidence decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep
[org_loc, org_conf] = deal(NaN(numel(sA), numel(sV), n_cue, numel(sigVs), n_rep));

fixP.sA = sA(a);
fixP.sV = sV(v);
fixP.sim_d = sim_d;
fixP.sigMotor = 1.36; % in deg, measured from first four participants
fixP.bi_nrep = n_rep;

[org_loc, org_conf] = sim_biloc_conf_5D(...
    aA, bA, sigA, sigV1, sigV2, muP, sigP, sigConf, pCommon, c1, c2, c3, lapse, fixP);

