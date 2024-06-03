% This model fits the unimodal localization and confidence responses in the
% first stage, and use the same \sigma_c, c1, c2, c3 to fit bimodal data.

clear; close all; rng('shuffle');

sub_slc = 12;
fitUniBi = 1;
if fitUniBi; ses_slc = 1:2; end % bimodal sessions

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'uniBiLoc_2stage')))

%% set cores

numCores = feature('numcores'); % number of cores locally
fprintf('Number of cores: %i  \n', numCores);
num_run = numCores - 1; % runs of fitting

%% part 1: fit unimodal loc and conf

% organize data
% condition (A,V1,V2) x loc (4) x rep
[data.uni_resp, data.uni_conf, ~, data.uniExpInfo] = org_data(sub_slc,[],'uniLoc');

% motor noise
data.sigMotor = get_point_sigM(sub_slc);

% general setting for all modelsmodel.num_run         = num_run;
model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run
model.x               = (-512:1:512) * deg_per_px;
model.sA              = data.uniExpInfo.audLevel;
model.sV              = model.sA;
model.num_rep         = data.uniExpInfo.nRep;

OPTIONS.TolMesh = 1e-5;
OPTIONS.Display = 'off';

% initiate
model.mode                  = 'initiate';
Val = nllUniLocConf([], model, data);

% optimize
model.mode                  = 'optimize';
NLL                         = NaN(1, model.num_runs);
estP                        = NaN(model.num_runs, Val.num_para);

parfor i                    = 1:model.num_runs

    disp(i);
    tempModel            = model;
    tempVal              = Val;
    tempFunc             = currModel;

    [estP(i,:),NLL(i)] = bads(@(p) tempFunc(p, model, data),...
        Val.init(i,:), Val.lb,...
        Val.ub, Val.plb, Val.pub, [], OPTIONS);

    disp(estP(i,:))

end














% sA(4) x sV(4) x post-cue(A, V) x reliability(high, low) x rep
[data.bi_resp, data.bi_conf, ~, data.biExpInfo] = org_data(sub_slc,ses_slc,'biLoc');