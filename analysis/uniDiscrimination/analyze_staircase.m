clear; clc; close all;

%% set up

sub = 'LL';
ses = 'A';
save_fig = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir, 'data','uniDiscrimination')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load data
load(sprintf('uniDis_sub-%s_ses-%s', sub, ses));
n_staircase = ExpInfo.n_staircase;

% replace nan trials
Resp.discrepancy(1,39:40) = -7.5;
Resp.discrepancy(2,39:40) = 7.5;
Resp.correct(1,39:40) = 1;
Resp.correct(2,39:40) = 1;

%% sort data

bin_width = 1;
edges =(min(ExpInfo.comparison_loc)-bin_width/2:bin_width:max(ExpInfo.comparison_loc)+bin_width/2) - mean(ExpInfo.standard_loc);

for ss = 1:n_staircase
    [sorted_discrepancy, sorted_idx] = sort(Resp.discrepancy(ss,:));
    org_r = Resp.correct(ss, sorted_idx);

    trial_counts = histcounts(sorted_discrepancy, edges); 
    correct_counts = histcounts(org_r, edges); 
end


%% plot the interleaved staircases


%% plot raw data for PMF

