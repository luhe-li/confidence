clear; close all; rng('shuffle');

%% set environment

useCluster                  = true;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

switch useCluster
    case true

        % See how many cores we have:
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores                    = maxNumCompThreads;
        end
        fprintf('Number of cores %i  \n', numCores);

        % Make sure Matlab does not exceed this
        maxNumCompThreads(numCores);

        hpc_job_number              = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end

        % set each job for a simulating model
        numJob                      = numel(hpc_job_number);
        fprintf('job number: %i \n', hpc_job_number);
        simd                       = hpc_job_number;

        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false

        % number of cores locally
        numCores = feature('numcores'); % number of cores locally
        fprintf('Number of cores: %i  \n', numCores);
        simd = 3;

end

%% set recovery parameters

% set experimental repetition per condition
num_rep = 30;
num_run = numCores - 1;

sampleGT = 1; % use a fixed set of GT or sample num_sample from GT range
if sampleGT
    num_sample = 100; checkFakeData = false;
else
    num_sample = 1; checkFakeData = true;
end

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV3_4ratings')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% experiment info

% fix parameters for the experiment info
speaker_span          = 65.5 * 2; % cm
sitting_dist          = 113; % cm
screen_width          = 170; % cm
screenX               = 1024; % pixel
screenXdeg            = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid            = screenX ./2;
num_speaker_int       = 15; % 15 intervals between 16 speakers
aud_level             = [6 8 9 11];
aud_VA                = -30:4:30;
deg_per_px            = screenXdeg / screenX;

sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA; % assume perceptually aligned
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of conditions

diff                 = zeros(length(sA), length(sV));
for ii                 = 1:length(sA)
    for j                 = 1:length(sV)
        diff(ii, j)           = sA(ii) - sV(j);
    end
end
all_diffs             = unique(abs(diff))';

%% model info

% localization parameter GT
if ~sampleGT
    %        aA,    bA, sigV1,  sigA,  sigV2, sigP,   pCC,  sigC
    GT = {[   1.1,   0.1,     1,   3,    4,    8,  0.6,   0.5],...% Heuristic
        [     1.1,   0.1,     1,   3,    4,    8,  0.6,   0.5],...% Suboptimal
        [     1.1,   0.1,     1,   3,    4,    8,  0.6,   0.5]}; % Optimal

else

    % TO-DO: values measured from fitting from pilot data
    mu_GT = [ 1.1,   0.1,     1,   3,    4,    8,  0.7,   0.5];
    se_GT = [0.21, 1.5, 0.19, 0.63, 0.8, 1.85, 0.1, 0.1];
    GT_samples = sampleGTfromGaussian(mu_GT, se_GT, num_sample);

end

ds_loc                = {'Model averaging'};
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

%% sample and fit for each simulating model

flnm = sprintf('RecoveryResults_simM%i_rep%i_sample%i_run%i.mat', simd, num_rep, num_sample, num_run);

[saveData, saveLocModel, saveConfModel] = deal(cell(num_model, num_model, num_sample));
[loc_CM, conf_CM] = deal(zeros(3));

for i_sample = 1:num_sample

    % specific GT for this sample
    if ~sampleGT
        i_gt = GT{i_sample, simd};
    else
        i_gt = GT_samples(i_sample,:);
    end

    % assign simulation parameters
    num_para              = length(i_gt);
    aA                    = i_gt(1);
    bA                    = i_gt(2);
    sigV1                 = i_gt(3);
    sigA                  = i_gt(4);
    sigV2                 = i_gt(5);
    sigVs                 = [sigV1, sigV2];
    sigP                  = i_gt(6);
    pCommon               = i_gt(7);
    sigC                  = i_gt(8);
    lapse                 = 0.02;
    muP                   = 0;

    %% set criterion based on other free parameters

    M_fixP.sA = sA;
    M_fixP.sV = sV;
    M_fixP.model_ind = simd;
    M_fixP.num_rep = num_rep;

    [lb, ub] = findConfRange(aA, bA, sigA, sigV1, sigV2, sigP, pCommon, M_fixP);
    intervals = linspace(lb, ub, 5);
    c1 = intervals(2);
    c2 = intervals(3);
    c3 = intervals(4);

    %% simulate fake data

    % num_s: stimulus location combination
    % 3 confidence decision strategies
    % 2 modalities(1 = aud, 2 = vis)
    % 2 visual reliabilities
    % num_rep
    [org_loc, org_conf] = deal(NaN(numel(sA), numel(sV), num_cue, numel(sigVs), num_rep));

    for a = 1:numel(sA)

        for v = 1:numel(sV)

            for r = 1:numel(sigVs)

                fixP.sA = sA(a);
                fixP.sV = sV(v);
                fixP.model_ind = simd;
                fixP.sigMotor = 1.36; % in deg, emperical motor noise averaged from first four participants
                fixP.num_rep = num_rep;

                [org_loc(a,v,:,r,:), org_conf(a,v,:,r,:)] = simAllModels(...
                    aA, bA, sigA, sigVs(r), muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP);

            end

        end

    end

    %% check fake data

    if checkFakeData

        %% check localization data

        uni_loc = zeros(size(org_loc));

        loc_a = repmat((sA*aA+bA)',[1,numel(sV)]);
        loc_v = repmat(sV,[numel(sA),1]);

        uni_loc(:,:,1,1:2,:) = repmat(loc_a, [1, 1, 1, 2, num_rep]);
        uni_loc(:,:,2,1:2,:) = repmat(loc_v, [1, 1, 1, 2, num_rep]);

        % loc at uni minus loc at bi
        ve =  mean(org_loc,5) - mean(uni_loc, 5);

        % diff x cue x reliability
        [ve_by_raw_diff, raw_diff] = org_by_raw_diffs_4D(ve, sA);

        % assume participants localized perfectly in the unisensory
        % condition
        figure; hold on
        t = tiledlayout(1, 2);
        title(t,sprintf('%s, rep: %i', ds_conf{simd}, num_rep))
        xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
        ylabel(t, 'Shift of localization');
        t.TileSpacing = 'compact';
        t.Padding = 'compact';

        for cue = 1:num_cue
            nexttile
            title(cue_label{cue})
            axis equal
            hold on

            for rel = 1: numel(sigVs)

                i_ve = squeeze(ve_by_raw_diff(:, cue, rel));
                plot(raw_diff, i_ve, 'Color',clt(rel+1,:))

            end
            xlim([min(raw_diff), max(raw_diff)])
            xticks(raw_diff)
            yline(0,'--')
            if cue == 1
                plot(raw_diff, raw_diff,'k--')
            else
                plot(raw_diff, -raw_diff,'k--')
            end
        end

        %% check confidence data

        % organize confidence data: {diff} cue x reliability x rep
        [conf_by_diff, all_diffs] = org_by_diffs(org_conf, sA);

        figure; hold on
        t = tiledlayout(2, 1);
        title(t,sprintf('%s, rep: %i', ds_conf{simd}, num_rep))
        xlabel(t, 'Absolute audiovisual discrepancy (deg)');
        ylabel(t, 'Confidence rating (1-4)');
        t.TileSpacing = 'compact';
        t.Padding = 'compact';

        for cue = 1:num_cue
            nexttile
            title(cue_label{cue})
            hold on
            for rel = 1: numel(sigVs)
                [p_conf, se_conf] = deal(NaN(1, numel(all_diffs)));
                for diff = 1:numel(all_diffs)
                    i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
                    p = sum(i_conf)/(numel(i_conf));
                    p_conf(diff) = p;
                    se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
                end
                plot(all_diffs, p_conf, 'Color',clt(rel+1,:));
                ylim([1, 4]) % rating range
            end
            xticks(all_diffs)
        end

    end

    %% fit the same data by all models

    fprintf('[%s] Start fitting simulating model-%i, sampled dataset-%i\n', mfilename, simd, i_sample);

    for fitd = 1:3

        % general setting for all models
        model.num_run         = num_run;
        model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run
        model.sA              = sA;
        model.sV              = model.sA;
        model.num_rep         = num_rep;
        model.num_SD          = 5;
        model.numBins_A       = 15;
        model.numBins_V       = 15;
        model.modality        = {'A','V'};
        model.strategy_loc    = 'MA';

        OPTIONS.TolMesh = 1e-5;
        OPTIONS.Display = 'off';
        %                 Turn limits on for formal fitting
        %                 OPTIONS.MaxIterations = 10000;
        %                 OPTIONS.MaxFunctionEvaluations = 10000;

        data.gt               = [i_gt, c1, c2-c1, c3-c2];
        data.bi_loc           = org_loc;
        data.bi_conf          = org_conf;
        data.sigMotor         = fixP.sigMotor;
        saveData{simd,fitd,i_sample}         = data;

        %% 1. first part: loc only

        % localization only
        currModel = str2func('nllLoc');

        % switch confidence strategies
        model.strategy_conf         = ds_conf{fitd};

        % initiate
        model.mode                  = 'initiate';
        Val = currModel([], model, data);

        % optimize
        model.mode                  = 'optimize';
        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

        parfor n              = 1:model.num_run

            tempModel             = model;
            tempVal               = Val;
            tempFunc              = currModel;

            [estP(n,:),NLL(n),~,~,~]    = bads(@(p) tempFunc(p, model, data),...
                Val.init(n,:), Val.lb,...
                Val.ub, Val.plb, Val.pub, [], OPTIONS);

        end

        % find the parameter with the least NLL
        [minNLL, best_idx]    = min(NLL);
        bestP                 = estP(best_idx, :);

        % save all fitting results
        saveLocModel{simd,fitd, i_sample}.paraInfo = Val;
        saveLocModel{simd,fitd,i_sample}.estP = estP;
        saveLocModel{simd,fitd,i_sample}.NLL = NLL;
        saveLocModel{simd,fitd,i_sample}.bestP = bestP;
        saveLocModel{simd,fitd,i_sample}.minNLL = minNLL;

        %% 2. second part: loc + conf

        % use best-fitting parameter to find out range for criteria
        M_fixP.sA = sA;
        M_fixP.sV = sV;
        M_fixP.model_ind = fitd;
        M_fixP.num_rep = num_rep;
        [lb, ub] = findConfRange(bestP(1), bestP(2), bestP(4), bestP(3), bestP(5), bestP(6), bestP(7), M_fixP);
        model.c_lb = lb;
        model.c_ub = ub;

        % use best-fitting parameter to set range for sigma_p
        mu_sig_p = mean(saveLocModel{simd,fitd,i_sample}.estP(:,6));
        sd_sig_p = std(saveLocModel{simd,fitd,i_sample}.estP(:,6));
        model.sig_p_lb = mu_sig_p - 3 * sd_sig_p;
        model.sig_p_ub = mu_sig_p + 3 * sd_sig_p;

        % loc + conf
        currModel = str2func('nllLocConf');

        % switch confidence strategies
        model.strategy_conf         = ds_conf{fitd};

        % initiate
        model.mode                  = 'initiate';
        Val = currModel([], model, data);

        % optimize
        model.mode                  = 'optimize';
        NLL                         = NaN(1, model.num_run);
        estP                        = NaN(model.num_run, Val.num_para);

        %                 p = [i_gt, c1, c2, c3];
        %                 test = currModel(p, model, data);

        parfor n              = 1:model.num_run

            tempModel             = model;
            tempVal               = Val;
            tempFunc              = currModel;

            [estP(n,:),NLL(n),~,~,~]    = bads(@(p) tempFunc(p, model, data),...
                Val.init(n,:), Val.lb, Val.ub, Val.plb, Val.pub, [], OPTIONS);

        end

        % find the parameter with the least NLL
        [minNLL, best_idx]    = min(NLL);
        bestP                 = estP(best_idx, :);

        % save all fitting results
        saveConfModel{simd,fitd,i_sample}.paraInfo = Val;
        saveConfModel{simd,fitd,i_sample}.estP = estP;
        saveConfModel{simd,fitd,i_sample}.NLL = NLL;
        saveConfModel{simd,fitd,i_sample}.bestP = bestP;
        saveConfModel{simd,fitd,i_sample}.minNLL = minNLL;

    end

    % find the best fitting model for this sampled dataset
    [~, best_fitd] = min(saveLocModel{simd,1,i_sample}.minNLL, saveLocModel{simd,2,i_sample}.minNLL, saveLocModel{simd,3,i_sample}.minNLL);
    loc_CM(simd, best_fitd) = loc_CM(simd, best_fitd) + 1;
    [~, best_fitd] = min(saveConfModel{simd,1,i_sample}.minNLL, saveConfModel{simd,2,i_sample}.minNLL, saveConfModel{simd,3,i_sample}.minNLL);
    conf_CM(simd, best_fitd) = conf_CM(simd, best_fitd) + 1;

    % save partial results
    save(fullfile(out_dir, flnm), 'saveLocModel','saveData','saveConfModel','loc_CM','conf_CM');

end

% save full results
fprintf('[%s] Model recovery done! Saving full results.\n', mfilename);
save(fullfile(out_dir, flnm), 'saveLocModel','saveData','saveConfModel','loc_CM','conf_CM');


%% utility

function [samples] = sampleGTfromGaussian(mean_values, sem_values, num_sample)

% Initialize the matrix to store the sampled values
samples = zeros(num_sample, 8);

% Parameters 'aA' and 'bA' are normally distributed
samples(:, 1) = normrnd(mean_values(1), sem_values(1), [num_sample, 1]);
samples(:, 2) = normrnd(mean_values(2), sem_values(2), [num_sample, 1]);

% Parameters '\sigma_{V1}', '\sigma_{A}', '\sigma_{V2}', '\sigma_{P}', '\sigma_{C}' are log-normally distributed
m = [mean_values(3), mean_values(4), mean_values(5), mean_values(6), mean_values(8)];
v = [sem_values(3:6), sem_values(8)].^2; % variance
log_mu = log((m.^2)./sqrt(v+m.^2));
log_sigma = sqrt(log(v./(m.^2)+1));
samples(:, 3) = lognrnd(log_mu(1), log_sigma(1), [num_sample, 1]);
samples(:, 4) = lognrnd(log_mu(2), log_sigma(2), [num_sample, 1]);
samples(:, 5) = lognrnd(log_mu(3), log_sigma(3), [num_sample, 1]);
samples(:, 6) = lognrnd(log_mu(4), log_sigma(4), [num_sample, 1]);
samples(:, 8) = lognrnd(log_mu(5), log_sigma(5), [num_sample, 1]);

% Parameter 'p_{common}' should be between 0 and 1, sampled from a truncated normal distribution
p_common_samples = normrnd(mean_values(7), sem_values(7), [num_sample, 1]);
p_common_samples(p_common_samples < 0) = 0;
p_common_samples(p_common_samples > 1) = 1;
samples(:, 7) = p_common_samples;

end