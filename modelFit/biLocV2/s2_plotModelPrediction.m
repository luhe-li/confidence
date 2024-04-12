
clear; clc; close all; rng('Shuffle');

sub = 4;
ses = 1:2;
models = {'Heuristic','Suboptimal','Optimal'};
m = 3;

flnm = sprintf('FitResults_%s_sub%i_ses%i-%i', models{m}, sub, min(ses), max(ses));

%% manage path
out_dir               = fullfile(pwd, mfilename);

% look for all the files, load the most updated one that match the selected name
files = dir(fullfile(pwd, [flnm '*.mat'])); % Lists all .mat files in the current folder
load(files(end).name);

% use the best-fitting parameter, simulate responses
p = model.bestP;

%% simulation

[loc, conf] = deal(NaN(num_s, num_cue, numel(sigVs), num_rep));
% num_s: stimulus location combination
% 3 confidence decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep

for j                 = 1:num_s
    for v                 = 1:numel(sigVs)
        [loc(j,:,v,:), conf(j,:,v,:)] = sim_loc_conf(num_rep, sAV(1,j), sAV(2,j),...
            aA, bA, sigA, sigVs(v), muP, sigP, pCommon, sigM, cA, cV, lapse, fixP, d);
    end
end

% organize predicted data


% plot

