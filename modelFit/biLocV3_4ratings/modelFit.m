clear; close all; rng('shuffle');

sub_slc = 12;
ses_slc = 1; % bimodal sessions

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))

%% set cores

numCores = feature('numcores'); % number of cores locally
fprintf('[%s] Number of cores: %i  \n', mfilename, numCores);
num_run = numCores - 1; % runs of fitting

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

    % sA(4) x sV(4) x post-cue(A, V) x reliability(high, low) x rep
    [data.bi_loc, data.bi_conf, ~, data.biExpInfo] = org_data(sub_slc,ses_slc,'biLoc');

    % motor noise
    data.sigMotor = get_point_sigM(sub_slc);

    % bi expt info for model
    model.bi_sA         = data.biExpInfo.speakerLocVA(data.biExpInfo.audIdx);
    model.bi_sV         = unique(data.biExpInfo.randVisVA);
    model.bi_nrep       = size(data.bi_loc, 5); % total number of trials per condition

    for d = 3:-1:1

        fprintf('[%s] Start fitting sub-%i, ses1-%i, M%i\n', mfilename, sub_slc, max(ses_slc), d);

        %% 1. first part: loc only

        % switch confidence strategies
        model.model_slc             = d;
        model.strategy_conf         = models{d};

        % initiate
        model.mode                  = 'initiate';
        Val = nllLoc([], model, data);

        % optimize
        model.mode                  = 'optimize';
        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

        parfor n              = 1:model.num_run

            tempModel             = model;
            tempVal               = Val;

            [estP(n,:),NLL(n),~,~,~]    = bads(@(p) nllLoc(p, model, data),...
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

        %% 2. second part: loc + conf

        % use best-fitting parameter to find out range for criteria
        M_fixP.sA = model.bi_sA;
        M_fixP.sV = model.bi_sV;
        M_fixP.model_ind = d;
        M_fixP.num_rep = model.bi_nrep;
        [lb, ub] = findConfRange(bestP(1), bestP(2), bestP(4), bestP(3), bestP(5), bestP(6), bestP(7), M_fixP);
        model.c_lb = lb;
        model.c_ub = ub;

        % use best-fitting parameter to set range for sigma_p
        mu_sig_p = mean(saveLocModel{d}.estP(:,6));
        sd_sig_p = std(saveLocModel{d}.estP(:,6));
        model.sig_p_lb = mu_sig_p - 3 * sd_sig_p;
        model.sig_p_ub = mu_sig_p + 3 * sd_sig_p;

        % initiate
        model.mode                  = 'initiate';
        Val = nllLocConf([], model, data);

        % optimize
        model.mode                  = 'optimize';
        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

        %                 p = [i_gt, c1, c2, c3];
        %                 test = currModel(p, model, data);

        parfor n              = 1:model.num_run

            tempModel             = model;
            tempVal               = Val;

            [estP(n,:),NLL(n),~,~,~]    = bads(@(p) nllLocConf(p, model, data),...
                Val.init(n,:), Val.lb, Val.ub, Val.plb, Val.pub, [], OPTIONS);

        end

        % find the parameter with the least NLL
        [minNLL, best_idx]    = min(NLL);
        bestP                 = estP(best_idx, :);

        % save all fitting results
        saveConfModel{d}.modelInfo = model;
        saveConfModel{d}.paraInfo = Val;
        saveConfModel{d}.estP = estP;
        saveConfModel{d}.NLL = NLL;
        saveConfModel{d}.bestP = bestP;
        saveConfModel{d}.minNLL = minNLL;

    end

    flnm = sprintf('fitResults_sub%i_ses%i-%i', sub, min(ses_slc), max(ses_slc));
    save(fullfile(out_dir, flnm), 'saveLocModel','data','saveConfModel','model');

end