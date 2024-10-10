clear; clc; close all;

%% set up

sub = 'OY';
save_fig = 0;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, mfilename);
addpath(genpath(fullfile(project_dir, 'util')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean

[org_uni, exp_info_uni] = org_resp(sub, {'A','V'}, 'uniLoc');
[org_bi, exp_info_bi] = org_resp(sub, 1, 'biLoc');
org = struct();
for i = 1:length(fieldnames(org_uni))
    org.(fields_uni{i}) = org_uni.(fields_uni{i});
end
for i = 1:length(fieldnames(org_bi))
    org.(fields_bi{i}) = org_bi.(fields_bi{i});
end


