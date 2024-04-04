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

%% start
flnm = sprintf('recoveryResults_rep%i_sample%i_run%i', num_rep, num_sample, num_runs);
