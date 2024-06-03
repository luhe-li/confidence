% This model fits the unimodal localization and confidence responses in the
% first stage, and use the same \sigma_c, c1, c2, c3 to fit bimodal data.

clear; close all; rng('shuffle');

sub_slc = 12;
fitUniBi = 1;
if fitUniBi; ses_slc = 1; end % bimodal sessions

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

%% organize data
% condition (A,V1,V2) x loc (4) x rep
[~, ~, ~, data.uniExpInfo, ~, ~, data.uni_loc, data.uni_conf] = org_data(sub_slc,[],'uniLoc');

% sA(4) x sV(4) x post-cue(A, V) x reliability(high, low) x rep
[data.bi_loc, data.bi_conf, ~, data.biExpInfo] = org_data(sub_slc,ses_slc,'biLoc');

% motor noise
data.sigMotor = get_point_sigM(sub_slc);

%% joint fit

flnm = sprintf('fitResults_sub%i_bises%i_run%i.mat', sub_slc, max(ses_slc), num_run);

% general setting for all models
model.num_run         = num_run;
model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run

% uni expt info for model
model.uni_sA          = data.uniExpInfo.speakerLocVA(data.uniExpInfo.audLevel);
model.uni_sV          = model.uni_sA;
model.uni_nrep        = data.uniExpInfo.nRep;

% bi expt info for model
model.bi_sA          = data.biExpInfo.speakerLocVA(data.biExpInfo.audIdx);
model.bi_sV          = unique(data.biExpInfo.randVisVA);
model.uni_nrep       = size(data.bi_loc, 5); % total number of trials per condition

% model setting for bimodal fitting
model.num_SD          = 5;
model.numBins_A       = 15;
model.numBins_V       = 15;
model.modality        = {'A','V'};
model.strategy_loc    = 'MA';
OPTIONS.TolMesh = 1e-5;
% OPTIONS.Display = 'off';

% confidence model setting
ds_conf               = {'Heuristic','Suboptimal','Optimal'};

% initiate
model.mode                  = 'initiate';
Val = nllUniBiLocConf([], model, data);

% optimize
model.mode                  = 'optimize';


for d = 3%:-1:1%1:3

    % switch confidence strategies
    model.strategy_conf         = ds_conf{d};

%     % test
%     p = [1:14]./14;
%     test = nllUniBiLocConf(p, model, data);

    NLL                         = NaN(1, model.num_run);
    estP                        = NaN(model.num_run, Val.num_para);

    parfor i                    = 1:model.num_run

        [estP(i,:),NLL(i)] = bads(@(p) nllUniBiLocConf(p, model, data),...
            Val.init(i,:), Val.lb,...
            Val.ub, Val.plb, Val.pub, [], OPTIONS);

    end

    % find the parameter with the least NLL
    [minNLL, best_idx]    = min(NLL);
    bestP                 = estP(best_idx, :);

    % save all fitting results
    saveModel{d}.paraInfo = Val;
    saveModel{d}.estP = estP;
    saveModel{d}.NLL = NLL;
    saveModel{d}.bestP = bestP;
    saveModel{d}.minNLL = minNLL;

    save(fullfile(out_dir, flnm), 'saveModel')

end
