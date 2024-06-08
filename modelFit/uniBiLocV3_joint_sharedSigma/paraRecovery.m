
clear; close all; rng('shuffle');

% sample participant
sub_slc = 15;
ses_slc = 1:3;

recompute = true;
sampleGT = 1; % use a fixed set of GT or sample num_sample from GT range

if sampleGT
    num_sample = 100; checkFakeData = false;
else
    num_sample = 1; checkFakeData = true;
end

%% set environment

useCluster                  = false;

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

        numJob                      = numel(hpc_job_number);
        fprintf('job number: %i \n', hpc_job_number);
        i_rep                       = hpc_job_number;

        % set experimental repetition per condition
        reps = round(logspace(log10(10), log10(100), 10));
        num_rep = reps(i_rep);

        % set run of fit by number of cores
        num_run = numCores - 1;

        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false
        numCores = feature('numcores'); % number of cores locally
        fprintf('Number of cores: %i  \n', numCores);
        num_rep = 100;
        num_run = numCores - 1;
end


%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))

%% experimental setup

% general setting for all models
model.num_run         = num_run;
model.num_sec         = num_run*2; % number of samples in the parameter space, must be larger than num_run
model.num_SD          = 5;
model.numBins_A       = 15;
model.numBins_V       = 15;
model.modality        = {'A','V'};
model.strategy_loc    = 'MA';
OPTIONS.TolMesh = 1e-5;

% info from data
[~, ~, ~, data.uniExpInfo, ~, ~, ~, ~, data.coefsA] = org_data(sub,[],'uniLoc');
[~, ~, ~, data.biExpInfo] = org_data(sub,1,'biLoc');

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

% model info
models               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(models);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

%% run

flnm = sprintf('recoveryResults_rep%i_sample%i_run%i.mat', num_rep, num_sample, num_run);

if ~exist(fullfile(out_dir, flnm),'file') || recompute

    [saveData, saveLocModel, saveModel, pred] = deal(cell(num_model, num_sample));

    for d = 3:-1:1

        %% set GT

        % load best-fitting parameters
        flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
        files = dir(fullfile(result_dir, flnm));
        F = load(files(end).name);

        % winning model?
        GTs(d,:) = F.saveModel{d}.bestP;

        if ~sampleGT
            GT_samples = fData.GTs(d,:);
        else
            mu_GT = GTs(d,:);
            se_GT = [0.21, 1.5, 0.19, 0.63, 0.8, 1.85, 0.1, 0.1];
            GT_samples = sampleGTfromGaussian(mu_GT, se_GT, num_sample);
        end

        %% simulate data

        model.mode = 'predict';
        model.model_slc             = d;
        model.strategy_conf         = models{d};

        for i_sample = 1:num_sample

            pred = nllUniBiLocConf(GT_samples(i_sample,:), model, data);
            fData{d,i_sample}.pred = pred;
            fData{d,i_sample}.GTs = GT_samples;

            %% fit by the same model

            fprintf('[%s] Start fitting sample-%i, rep-%i, M%i\n', mfilename, i_sample, num_rep, d);

            % initiate
            model.mode                  = 'initiate';
            Val = nllUniBiLocConf([], model, data);

            % optimize
            model.mode                  = 'optimize';

            NLL                         = NaN(1, model.num_run);
            estP                        = NaN(model.num_run, Val.num_para);

            parfor i_run                    = 1:model.num_run

                [estP(i_run,:),NLL(i_run)] = bads(@(p) nllUniBiLocConf(p, model, data),...
                    Val.init(i_run,:), Val.lb,...
                    Val.ub, Val.plb, Val.pub, [], OPTIONS);

            end

            % find the parameter with the least NLL
            [minNLL, best_idx]    = min(NLL);
            bestP                 = estP(best_idx, :);

            % save all fitting results
            saveModel{d,i_sample}.paraInfo = Val;
            saveModel{d,i_sample}.estP = estP;
            saveModel{d,i_sample}.NLL = NLL;
            saveModel{d,i_sample}.bestP = bestP;
            saveModel{d,i_sample}.minNLL = minNLL;

        end

        save(fullfile(out_dir, flnm), 'saveModel','model','saveConfModel')
    end



end


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