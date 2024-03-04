clear; close all;

%% key model recovery parameters

num_rep               = 12;
num_runs              = 10;
num_sample            = 10;

%% Manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, 'v1_parameter_recovery');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV1')))

%% Experimental info

% fix parameters for the experiment info
speaker_span          = 65.5 * 2; % cm
sitting_dist          = 113; % cm
screen_width          = 170; % cm
screenX               = 1024; % pixel
screenXdeg            = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid            = screenX ./2;
num_speaker_int       = 15; % 15 intervals between 16 speakers
cm_per_aud_ind        = speaker_span / num_speaker_int;
pixel_per_cm          = screenX / screen_width;
aud_level             = [6 8 9 11];
aud_VA                = -30:4:30;
deg_per_px            = screenXdeg / screenX;

fixP.screenX          = screenXdeg;
fixP.x                = -screenXdeg /2 : deg_per_px : screenXdeg/2; % the screen deg space
center_x              = 0;

sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA;
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of conditions

aud_locs              = sA;
vis_locs              = sV;
diffs                 = zeros(length(aud_locs), length(vis_locs));
for i                 = 1:length(aud_locs)
    for j                 = 1:length(vis_locs)
        diffs(i, j)           = aud_locs(i) - vis_locs(j);
    end
end
disc_locs             = unique(abs(diffs));

%% Free parameters

GT.aA                 = [   1,   1.5]; % scale
GT.bA                 = [-0.5,   0.5]; % intercept in degree
GT.sigA               = [ 0.5,     3]; % degree
GT.sigV1              = [ 0.1,     1]; % degree
GT.sigV2              = [ 0.5,     3]; % degree
GT.sigP               = [  10,    15]; % degree
GT.pC1                = [ 0.3,   0.7]; % weight
GT.c                  = [ 0.4,   0.6]; % weight

fn                    = fieldnames(GT);
num_para              = numel(fn);
[lb, ub]              = deal(NaN(1,  num_para));
for k                 = 1:numel(fn)
    lb(:,k)               = GT.(fn{k})(1);
    ub(:,k)               = GT.(fn{k})(2);
end

% randomly sample from the range
samples               = lb + (ub - lb) .* rand(num_sample, num_para);

% simulated model info
ds_loc                = {'Model averaging'};
ds_conf               = {'Optimal','Suboptimal','Heuristic'};
num_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

%% Simulate data

% preallocate
saveModel             = cell(1, num_model); % for each simulated moedl
for d                 = 1:num_model
    saveModel{d}.estP     = NaN(num_model, num_sample, num_runs, num_para); % first dimension is the fitting model
    saveModel{d}.NLL      = NaN(num_model, num_sample, num_runs);
    saveModel{d}.bestP    = NaN(num_model, num_sample, num_para);
    saveModel{d}.minNLL   = NaN(num_model, num_sample, 1);
end
pred                  = cell(num_sample, num_model, num_model);
CM                    = zeros(num_model);

for i                 = 1:num_sample

    disp(i)
    aA                    = samples(i, 1);
    bA                    = samples(i, 2);
    sigA                  = samples(i, 3);
    sigVs                 = samples(i, 4:5);
    sigP                  = samples(i, 6);
    pCommon               = samples(i, 7);
    criterion             = samples(i, 8);

    [loc, conf, variance] = deal(NaN(num_s,num_model, num_cue, numel(sigVs), num_rep));
    % num_s: stimulus location combination
    % 3 confidence decision strategies
    % 2 modalities(1 = aud, 2 = vis)
    % 2 visual reliabilities
    % num_rep

    for j                 = 1:num_s
        for v                 = 1:numel(sigVs)
            [loc(j,:,:,v,:), conf(j,:,:,v,:), variance(j,:,:,v,:)] = sim_loc_pconf(num_rep, sAV(1,j),...
                sAV(2,j), aA, bA, sigA, sigVs(v), sigP, pCommon, criterion, fixP);
        end
    end

    %% Organize data

    [org_resp, org_conf]  = deal(NaN(numel(ds_loc),numel(sA), numel(sV), num_cue, numel(sigVs), num_rep));

    for d                 = 1:num_model

        d_resp                = squeeze(loc(:,d,:,:,:));
        org_resp(d,:,:,:,:,:) = reshape(d_resp, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

        d_conf                = squeeze(conf(:,d,:,:,:));
        org_conf(d,:,:,:,:,:) = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

    end

    %% Model fitting

    % general setting for all models
    model.num_runs        = num_runs;
    model.x               = (-512:1:512) * deg_per_px;
    model.sA              = sA;
    model.sV              = model.sA;
    model.num_rep         = num_rep;

    % set OPTIONS to tell bads that my objective function is noisy
    OPTIONS.UncertaintyHandling                            = 1;


    for d                 = 1:num_model

        data.org_resp         = squeeze(org_resp(d,:,:,:,:,:));
        data.org_conf         = squeeze(org_conf(d,:,:,:,:,:));
        data.sigM             = 1.36; % emperical motor noise averaged from first four participants

        % fit by all models
        for k                = 1:num_model

            currModel             = str2func(['nll', ds_conf{k}]);

            model.mode            = 'initiate';
            Val                   = currModel([], model, data);

            model.mode            = 'optimize';
            NLL                   = NaN(1, model.num_runs);
            estP                  = NaN(model.num_runs, Val.num_para);

            parfor n              = 1:model.num_runs

                tempModel             = model;
                tempVal               = Val;
                tempFunc              = currModel;

                [estP(n,:),NLL(n)]    = bads(@(p) tempFunc(p, model, data),...
                    Val.init(n,:), Val.lb,...
                    Val.ub, Val.plb, Val.pub, [], OPTIONS);

                disp(estP(n,:))

            end

            % find the parameter with the least NLL
            [minNLL, best_idx]    = min(NLL);
            bestP                 = estP(best_idx, :);

            % save all fitting results
            saveModel{d}.estP(k,i,:,:) = estP;
            saveModel{d}.NLL(k,i,:) = NLL;
            saveModel{d}.bestP(k,i,:) = bestP;
            saveModel{d}.minNLL(k,i)= minNLL;

            % predict using the best-fitting parameter
            model.mode            = 'predict';
            tmpPred               = currModel(bestP, model, data);
            tmpPred.bestP         = bestP;
            pred{i, d, k}            = tmpPred;

        end

    end

end

save('recoveryResults')

%% Plot parameters (predicted vs. ground-truth)

for dd                = 1:num_model

    figure;
    t                     = tiledlayout(2, 4);
    title(t, sprintf('%s, rep: %i', ds_conf{dd}, num_rep));

    % Loop through each parameter
    for jj                = 1:8

        nexttile;

        % Scatter plot of the i-th predicted parameters against its ground-truth
        scatter(samples(jj,:), saveModel{d}.bestP(dd,:,jj), 'filled');
        xlabel('Ground Truth');
        ylabel('Predicted');
        axis square; % Make the plot square

        % Calculate the Pearson correlation coefficient and p-value
        [R, P]                = corrcoef(ground_truth(:, jj), saveModel{d}.bestP(dd,:,jj));

        % Extract the correlation coefficient and p-value
        r                     = R(1,2);
        p_value               = P(1,2);

        % Label the r and p in the title
        title(sprintf('Param %d: r= %.2f, p=%.3f', i, r, p_value));
    end

end

%% confusion matrix

CM = zeros(num_model);

for d = 1:num_model

    nll_of_sim_model = saveModel{d}.minNLL;

    % Step 1: Find the indices (rows) of the minimum values in each column
    [~, minIndices] = min(nll_of_sim_model);

    % Step 2: Calculate the percentage of each k element being the min
    minCounts = zeros(num_model, 1); 

    for idx = 1:num_sample
        minCounts(minIndices(idx)) = minCounts(minIndices(idx)) + 1;
    end

    minPercentages = (minCounts / num_sample)'; 
    CM(d,:)        = minPercentages;

end
