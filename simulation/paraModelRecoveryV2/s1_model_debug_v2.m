% - rewritten simulation codes and fitting codes
% - add data analysis checking (localization + confidence)

clear; close all;

num_rep = 1000; % rep of trials
num_run = 10; % run of fitting
num_sample = 1; % number of fake data samples
check_fake_data = 1;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV2')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;
%     251, 154, 153]./255; % light red


%% experiment info
speaker_span          = 65.5 * 2; % cm
sitting_dist          = 113; % cm
screen_width          = 170; % cm
screenX               = 1024; % pixel
screenXdeg            = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid            = screenX ./2;
num_speaker_int       = 15; % 15 intervals between 16 speakers
cm_per_aud_ind        = speaker_span / num_speaker_int;
pixel_per_cm          = screenX / screen_width;
aud_level             = [6 8 9 11];
aud_VA                = -30:4:30;
deg_per_px            = screenXdeg / screenX;

fixP.screenX          = screenXdeg;
fixP.x                = -screenXdeg /2 : deg_per_px : screenXdeg/2; % the screen deg space

sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA;
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of conditions

aud_locs              = sA;
vis_locs              = sV;
diffs                 = zeros(length(aud_locs), length(vis_locs));
for i                 = 1:length(aud_locs)
    for j                 = 1:length(vis_locs)
        diffs(i, j)           = aud_locs(i) - vis_locs(j);
    end
end
diff_locs             = unique(abs(diffs))';

%% start
flnm = sprintf('recoveryResults_rep%i_sample%i_run%i', num_rep, num_sample, num_run);

% choose a reasonble set of parameter set. See variable name below.
%    aA, bA, sigV1, dsigA, dsigV2, sigP,  pCC, sigM, cA, cV
GT = {[1,  0,  0.7,   1.2,    1.5,   15, 0.57,  0.5, 50, 100],...% Heuristic
      [1,  0,  0.7,   1.2,    1.5,   15, 0.57,  0.5, 50, 100],...% Suboptimal
      [1,  0,  0.7,   1.2,    1.5,   15, 0.57,  0.5, 50, 70]}; % Optimal

% simulated model info
ds_loc                = {'Model averaging'};
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

for d = 3:-1:1
    
    aA                    = GT(d,1);
    bA                    = GT(d,2);
    sigV1                 = GT(d,3);
    sigA                  = GT(d,4) + sigV1;
    sigV2                 = GT(d,5) + sigV1;
    sigVs                 = [sigV1, sigV2];
    sigP                  = GT(d,6);
    pCommon               = GT(d,7);
    sigM                  = GT(d,8);
    cA                    = GT(d,9);
    cV                    = GT(d,10);
    lapse                 = 0.02;
    muP                   = 0;
    
    
    [loc, conf] = deal(NaN(num_s, num_cue, numel(sigVs), num_rep));
    % num_s: stimulus location combination
    % 3 confidence decision strategies
    % 2 modalities(1 = aud, 2 = vis)
    % 2 visual reliabilities
    % num_rep
    
    for j                 = 1:num_s
        for v                 = 1:numel(sigVs)
            [loc(j,:,v,:), conf(j,:,v,:)] = sim_loc_conf_v2(num_rep, sAV(1,j), sAV(2,j),...
                aA, bA, sigA, sigVs(v), muP, sigP, pCommon, sigM, cA, cV, lapse,  fixP, d);
        end
    end
    
end
