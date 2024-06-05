% This model fits the unimodal localization and confidence responses in the
% first stage, and use the same \sigma_c, c1, c2, c3 to fit bimodal data.

clear; close all; rng('shuffle');

sub_slc = 15;
ses_slc = 1:3; % bimodal sessions

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'uniBiLoc_joint')))

%% set cores

numCores = feature('numcores'); % number of cores locally
fprintf('Number of cores: %i  \n', numCores);
num_run = numCores - 1; % runs of fitting

%% joint fit

% general setting for all models
model.num_run         = num_run;
model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run
model.num_SD          = 5;
model.numBins_A       = 15;
model.numBins_V       = 15;
model.modality        = {'A','V'};
model.strategy_loc    = 'MA';
OPTIONS.TolMesh = 1e-5;
% OPTIONS.Display = 'off';

models               = {'Heuristic','Suboptimal','Optimal'};

for sub = sub_slc

    %% organize data

    % condition (A,V1,V2) x loc (4) x rep
    [data.org_uni_loc, data.org_uni_conf, ~, data.uniExpInfo, ~, ~, data.uni_loc, data.uni_conf, data.coefsA] = org_data(sub,[],'uniLoc');

    % sA(4) x sV(4) x post-cue(A, V) x reliability(high, low) x rep
    [data.bi_loc, data.bi_conf, ~, data.biExpInfo] = org_data(sub,ses_slc,'biLoc');

    % motor noise
    data.sigMotor = get_point_sigM(sub);

    % uni expt info for model
    model.uni_sA          = unique(data.uniExpInfo.randAudVA);
    model.uni_sV          = unique(data.uniExpInfo.randVisVA);
    model.uni_nrep        = data.uniExpInfo.nRep;

    % bi expt info for model
    model.bi_sA          = unique(data.biExpInfo.randAudVA);
    model.bi_sV          = unique(data.biExpInfo.randVisVA);
    model.bi_nrep        = size(data.bi_loc, 5); % total number of trials per condition

    for d = 3:-1:1

        fprintf('[%s] Start fitting sub-%i, ses1-%i, M%i\n', mfilename, sub_slc, max(ses_slc), d);

        % switch confidence strategies
        model.strategy_conf         = models{d};
        % initiate
        model.mode                  = 'initiate';
        Val = nllUniBiLocConf([], model, data);

        % optimize
        model.mode                  = 'optimize';

%             % test
%             p = [ 1.7090   12.0101    5.6956    9.0483    0.4485    1.0075    1.7248    1.4398    1.3167    0.0600    5.0000];
%             test = nllUniBiLocConf(p, model, data);

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

        % save all fitting result
        saveModel{d}.paraInfo = Val;
        saveModel{d}.estP = estP;
        saveModel{d}.NLL = NLL;
        saveModel{d}.bestP = bestP;
        saveModel{d}.minNLL = minNLL;

    end

    flnm = sprintf('fitResults_sub%i_ses%i-%i', sub, min(ses_slc), max(ses_slc));
    save(fullfile(out_dir, flnm), 'saveModel','data','model');

end