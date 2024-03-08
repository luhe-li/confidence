
clear; clc; close all; rng('Shuffle');

%% set environment

useCluster                  = 0;
sub = 2;
ses = 1:2;

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

% set mode
model.mode                  = 'optimize';
model.num_runs              = 5;

% set fixed & set-up parameters
deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
model.x                     = (-512:1:512) * deg_per_px;
model.sA                    = ExpInfo.speakerLocVA(ExpInfo.audIdx);
model.sV                    = model.sA;

model.paraID                = {'aA','bA','\sigma_{A}','\sigma_{V1}','\sigma_{V2}','\sigma_{P}','p_{common}','criterion'};
model.num_para              = length(model.paraID);
model.num_rep               = size(data.org_resp, 5);

% hard bounds, the range for LB, UB, larger than soft bounds
paraH.aA                    = [ 0.5,     2]; % degree
paraH.bA                    = [  -8,     8]; % degree
paraH.sigA                  = [1e-2,    10]; % degree
paraH.sigV1                 = [1e-2,    10]; % degree
paraH.sigV2                 = [1e-2,    10]; % degree
paraH.sigP                  = [   1,    30]; % degree
paraH.pC1                   = [1e-4,1-1e-4]; % weight
paraH.c                     = [1e-4,1-1e-4]; % weight

% soft bounds, the range for PLB, PUB
paraS.aA                    = [ 0.8,   1.2]; % degree
paraS.bA                    = [  -4,     4]; % degree
paraS.sigA                  = [   4,     6]; % degree
paraS.sigV1                 = [   2,     6]; % degree
paraS.sigV2                 = [   4,     6]; % degree
paraS.sigP                  = [   1,    20]; % degree
paraS.pC1                   = [ 0.4,   0.6]; % weight
paraS.c                     = [ 0.2,   0.5]; % weight

% reorganize parameter bounds to feed to bads
fn                          = fieldnames(paraH);
for k                       = 1:numel(fn)
    model.lb(:,k)               = paraH.(fn{k})(1);
    model.ub(:,k)               = paraH.(fn{k})(2);
    model.plb(:,k)              = paraS.(fn{k})(1);
    model.pub(:,k)              = paraS.(fn{k})(2);
end
model.paraS                 = paraS; model.paraH = paraH;

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 0;

%% model fitting

% define function calculating nll
funcNLL                     = @(p) nllOptimal(p(1), p(2), p(3), p(4), p(5), p(6),p(7), p(8),...
    model, data);

% get grid initializations
numSections                 = model.num_runs;
model.init                  = getInit(model.lb, model.ub, numSections, model.num_runs);

%initialize matrices that store negative log likelihood and best-fit paramters
NLL                      = NaN(1, model.num_runs);
estP                  = NaN(model.num_runs, length(model.lb));

% test with a set of parameters if needed
p = [ 1.2, 4, 1.6992    0.0122    0.0342   29.9916    0.2557    0.4714];
testnll     = nllOptimal(p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), model, data);

for i                    = 1:model.num_runs
    disp(i);
    try
        tempModel                   = model;
        [estP(i,:),NLL(i)] = bads(funcNLL, tempModel.init(i,:), tempModel.lb,...
            tempModel.ub, tempModel.plb, tempModel.pub, [], OPTIONS);
        disp(estP(i,:))
    catch
        disp('Invalid NLL.')
        continue;
    end
end

model.estP            = estP;
model.NLL             = NLL;
[model.minNLL, best_idx] = min(NLL);
bestP = estP(best_idx, :);
model.bestP = bestP;
pred.bestP = bestP;

%% predict using the best-fitting parameter

p = bestP;
pred = nllOptimal(p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), model, data);

% save the data for each participant
flnm                        = sprintf('sub%d_ses%i-%i', sub, min(ses), max(ses));
save(['ModelFitResults_' flnm datestr(datetime('now'))],'data','model','pred')

