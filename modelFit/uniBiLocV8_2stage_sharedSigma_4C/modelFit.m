% This model fits the unimodal and bimodal localization data first, and
% freeze the parameters to fit confidence rating.

% Not the ideal way, but may help clarify necessary confidence model
% parameters

clear; close all; rng('shuffle');

sub_slc = 13;
ses_slc = 1:3; % bimodal sessions

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))

%% set cores

numCores = feature('numcores'); % number of cores locally
fprintf('Number of cores: %i  \n', numCores);
num_run = numCores - 5; % runs of fitting

%% model set up

% general setting for all models
model.num_run         = num_run;
model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run
model.num_SD          = 5;
model.numBins_A       = 15;
model.numBins_V       = 15;
model.modality        = {'A','V'};
model.strategy_loc    = 'MA';
OPTIONS.TolMesh = 1e-5;
OPTIONS.Display = 'off';

% loop each model
models               = {'Heuristic','Suboptimal','Optimal'};

for sub = sub_slc

    %% organize data

    % condition (A,V1,V2) x loc (4) x rep
    [data.org_uni_loc, data.org_uni_conf, ~, data.uniExpInfo, ~, ~, data.uni_loc, data.uni_conf] = org_data(sub,[],'uniLoc');

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

        fprintf('[%s] Start fitting localization data, sub-%i, ses1-%i, M%i\n', mfilename, sub_slc, max(ses_slc), d);

        % switch confidence strategies
        model.strategy_conf         = models{d};

        % initiate
        model.mode                  = 'initiate';
        Val = nllUniBiLoc([], model, data);

        % optimize
        model.mode                  = 'optimize';

        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

       
        parfor n              = 1:model.num_run

            tempModel             = model;
            tempVal               = Val;

            [estP(n,:),NLL(n),~,~,~]    = bads(@(p) nllUniBiLoc(p, model, data),...
                Val.init(n,:), Val.lb,...
                Val.ub, Val.plb, Val.pub, [], OPTIONS);

        end

        % find the parameter with the least NLL
        [minNLL, best_idx]    = min(NLL);
        bestP                 = estP(best_idx, :);

        % save all fitting results
        saveLocModel{d}.paraInfo = Val;
        saveLocModel{d}.estP = estP;
        saveLocModel{d}.NLL = NLL;
        saveLocModel{d}.bestP = bestP;
        saveLocModel{d}.minNLL = minNLL;

        %% 2. second part: conf

        fprintf('[%s] Start fitting confidence data, sub-%i, ses1-%i, M%i\n', mfilename, sub_slc, max(ses_slc), d);

        % use best param from loc fits
        model.locP = bestP;

        % initiate
        model.mode                  = 'initiate';
        Val = nllUniBiConf([], model, data);

        % optimize
        model.mode                  = 'optimize';

        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

        test = nllUniBiConf(Val.init(1,:), model, data);

        parfor i                    = 1:model.num_run

            [estP(i,:),NLL(i)] = bads(@(p) nllUniBiConf(p, model, data),...
                Val.init(i,:), Val.lb,...
                Val.ub, Val.plb, Val.pub, [], OPTIONS);

        end

        % find the parameter with the least NLL
        [minNLL, best_idx]    = min(NLL);
        bestP                 = estP(best_idx, :);

        % save all fitting result
        saveConfModel{d}.paraInfo = Val;
        saveConfModel{d}.estP = estP;
        saveConfModel{d}.NLL = NLL;
        saveConfModel{d}.bestP = bestP;
        saveConfModel{d}.minNLL = minNLL;

    end

    flnm = sprintf('fitResults_sub%i_ses%i-%i', sub, min(ses_slc), max(ses_slc));
    save(fullfile(out_dir, flnm), 'saveLocModel','data','saveConfModel','model');

end
