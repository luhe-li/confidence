
clear; clc; close all; rng('Shuffle');

%% set environment

useCluster                  = 0;
sub = 2;
ses = 1:2;
models = {'Optimal',...
    'Suboptimal',...
    'Heuristic'};

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
        %         sub                         = hpc_job_number;

        if isempty(gcp('nocreate'))
            parpool(numCores-1);
        end

    case false
        numCores                    = 8; % number of cores locally
        fprintf('Number of cores: %i  \n', numCores);
end

%% manage paths

currentDir                  = pwd;
[projectDir, ~]             = fileparts(fileparts(currentDir));
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'bads')));
addpath(genpath(fullfile(projectDir, 'func')));
addpath(genpath(fullfile(projectDir, 'exptCode/biloc/')));

%% organize data

[data.org_resp, data.org_conf, ~, ExpInfo, ~, ScreenInfo] = org_data(sub,ses,'biLoc');
data.sigM = get_point_sigM(sub);
load(fullfile(projectDir, ['exptCode/biloc/' sprintf('AVbias_sub%i', sub) '.mat']));
data.coefsA = squeeze(Transfer.degCoeff(1, :));

%% define model

% set fixed & set-up parameters
model.num_runs              = 7;
deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
model.x                     = (-512:1:512) * deg_per_px;
model.sA                    = ExpInfo.speakerLocVA(ExpInfo.audIdx);
model.sV                    = model.sA;
model.num_rep               = size(data.org_resp, 5);

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 1;

%% model fitting

for m = 1:numel(models)

    currModel = str2func(['nll', models{m}]);

    model.mode                  = 'initiate';
    Val = currModel([], model, data);

    model.mode                  = 'optimize';
    NLL                         = NaN(1, model.num_runs);
    estP                        = NaN(model.num_runs, Val.num_para);

    for i                    = 1:model.num_runs

        disp(i);
        tempModel            = model;
        tempVal              = Val;
        tempFunc             = currModel;

        [estP(i,:),NLL(i)] = bads(@(p) tempFunc(p, model, data),...
            Val.init(i,:), Val.lb,...
            Val.ub, Val.plb, Val.pub, [], OPTIONS);

        disp(estP(i,:))

    end

    model.estP            = estP;
    model.NLL             = NLL;

    % find the parameter with the least NLL
    [model.minNLL, best_idx] = min(NLL);
    bestP = estP(best_idx, :);
    model.bestP = bestP;

    % predict using the best-fitting parameter
    model.mode                  = 'predict';
    pred = currModel(bestP, model, data);
    pred.bestP = bestP;

    % save the data for each participant
    flnm        = sprintf('sub%d_ses%i-%i', sub, min(ses), max(ses));
    save(['FitResults_' models{m} '_' flnm datestr(datetime('now'))],'data','model','pred')

end
