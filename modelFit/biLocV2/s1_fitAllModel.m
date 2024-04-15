
clear; clc; close all; rng('Shuffle');

%% set environment

useCluster = 0;
sub = 4;
ses = 1:2;
models = {'Heuristic','Suboptimal','Optimal'};

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
data.sigMotor = get_point_sigM(sub);
% load(fullfile(projectDir, ['exptCode/biloc/' sprintf('AVbias_sub%i', sub) '.mat']));
% data.coefsA = squeeze(Transfer.degCoeff(1, :));

%% define model

% set fixed & set-up parameters
model.num_runs              = 3;
model.num_sec               = model.num_runs*2;
deg_per_px                  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
model.x                     = (-512:1:512) * deg_per_px;
model.sA                    = round(ExpInfo.speakerLocVA(ExpInfo.audIdx));
model.sV                    = model.sA;
model.num_SD                = 5;
model.numBins_A             = 10;
model.numBins_V             = 15;
model.modality              = {'A','V'};
model.num_rep               = size(data.org_resp, 5);
model.strategy_loc          = 'MA';

% set OPTIONS to tell bads that my objective function is noisy
OPTIONS.UncertaintyHandling = 0;
OPTIONS.TolMesh = 1e-4;

%% check data

    uni_loc = nan(size(data.org_resp));

    loc_a = repmat(model.sA',[1,numel(model.sV)]);
    loc_v = repmat(model.sV,[numel(model.sA),1]);

    uni_loc(:,:,1,:,:) = repmat(loc_a, [1, 1, 1, 2, model.num_rep]);
    uni_loc(:,:,2,:,:) = repmat(loc_v, [1, 1, 1, 2, model.num_rep]);

    % loc at uni minus loc at bi
    ve =  mean(data.org_resp,5) - mean(uni_loc, 5);

    % diff x cue x reliability
    [ve_by_raw_diff, all_raw_diffs] = org_by_raw_diffs_4D(ve, model.sA);

    % assume participants localized perfectly in the unisensory
    % condition
    figure; hold on
    t = tiledlayout(2, 1);
%     title(t,sprintf('M%s, rep: %i', model_slc, num_rep))
    xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
    ylabel(t, 'Shift of localization');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for cue = 1:2
        nexttile
        axis equal
        hold on

        for rel = 1: 2

            i_ve = squeeze(ve_by_raw_diff(:, cue, rel));
            plot(all_raw_diffs, i_ve)

        end
        xticks(all_raw_diffs)

        yline(0,'--')
        if cue == 1
            plot(all_raw_diffs, all_raw_diffs,'k--')
        else
            plot(all_raw_diffs, -all_raw_diffs,'k--')
        end
    end


%% model fitting

for m = numel(models):-1:1

    currModel = str2func('nllBimodal');

    % switch confidence strategies
    model.strategy_conf         = models{m};

    % initiate
    model.mode                  = 'initiate';
    Val = currModel([], model, data);

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

    model.estP            = estP;
    model.NLL             = NLL;

    % find the parameter with the least NLL
    [model.minNLL, best_idx] = min(NLL);
    bestP = estP(best_idx, :);
    model.bestP = bestP;

    % predict using the best-fitting parameter
    model.mode                  = 'predict';
    pred = currModel(bestP, model, data); % two structs for two reliabilities

    % save the data for each participant
    flnm        = sprintf('sub%d_ses%i-%i', sub, min(ses), max(ses));
    currentDateTimeStr = datestr(datetime('now'));
    modifiedDateTimeStr = strrep(currentDateTimeStr, ':', '-');
    save(['FitResults_' models{m} '_' flnm '_' modifiedDateTimeStr],'data','model','pred')

end
