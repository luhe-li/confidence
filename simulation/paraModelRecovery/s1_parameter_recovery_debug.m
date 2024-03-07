% simulate one set of parameter
% fit by v2: fitAllModels (know inside)

clear; close all;

%% key model recovery parameters

num_rep               = 12;
num_runs              = 7;
num_sample            = 1;
checkFakeData         = 1;

%% Manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, 's1_parameter_recovery');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV2')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;
%     251, 154, 153]./255; % light red

%% start
flnm = sprintf('recoveryResults_rep%i_sample%i_run%i', num_rep, num_sample, num_runs);

if exist(flnm,'file')

    load(flnm);

else

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
    diff_locs             = unique(abs(diffs))';

    %% Free parameters

    % choose a reasonble set of parameter set. See variable name below.

    GT = [zeros(1,9);...
        1, 0,   1,   4,10, 1e4, 0.01, 13.8, 13;...%Suboptimal
        1, 0, 0.7, 1.2, 1, 1e4, 0.57, 1.5, 0.55]; % Optimal
    num_para         = size(GT, 2);

    % simulated model info
    ds_loc                = {'Model averaging'};
    ds_conf               = {'Heuristic','Suboptimal','Optimal'};
    num_model             = numel(ds_conf);
    cue_label             = {'Post-cue: A','Post-cue: V'};
    num_cue               = numel(cue_label);
    rel_label             = {'High reliability','Low reliability'};

    %% Simulate data

    for d = 3:-1:1

        aA                    = GT(d,1);
        bA                    = GT(d,2);
        sigV1                 = GT(d,3);
        sigA                  = GT(d,4) + sigV1;
        sigV2                 = GT(d,5) + sigV1;
        sigVs                 = [sigV1, sigV2];
        sigP                  = GT(d,6);
        pCommon               = GT(d,7);
        cA                    = GT(d,8);
        cV                    = GT(d,9);
        lapse = 0.05;
        muP = 0;

        [loc, conf] = deal(NaN(num_s, num_cue, numel(sigVs), num_rep));
        % num_s: stimulus location combination
        % 3 confidence decision strategies
        % 2 modalities(1 = aud, 2 = vis)
        % 2 visual reliabilities
        % num_rep

        for j                 = 1:num_s
            for v                 = 1:numel(sigVs)
                [loc(j,:,v,:), conf(j,:,v,:)] = sim_loc_conf_2criteria(pCommon,...
                    num_rep, sAV(1,j), sAV(2,j), aA, bA, sigA, sigVs(v), muP, sigP, fixP, [cA, cV], lapse, d);
            end
        end

        %% Organize data

        org_resp              = reshape(loc, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);
        org_conf              = reshape(conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

        %% check fake data
        if checkFakeData

            %{diff} cue x reliability x rep
            [conf_by_diff, all_diffs] = org_by_diffs(org_conf, sA);

            figure; hold on
            t = tiledlayout(2, 1);
            title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
            xlabel(t, 'Audiovisual discrepancy (deg)');
            ylabel(t, 'Proportion of confidence report');
            t.TileSpacing = 'compact';
            t.Padding = 'compact';

            for cue = 1:num_cue
                nexttile
                title(cue_label{cue})
                hold on
                for rel = 1: numel(sigVs)
                    [p_conf, se_conf] = deal(NaN(1, numel(diff_locs)));
                    for diff = 1:numel(diff_locs)
                        i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
                        p = sum(i_conf)/numel(i_conf);
                        p_conf(diff) = p;
                        se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
                    end
                    plot(diff_locs, p_conf, 'Color',clt(rel+1,:));
                    ylim([-0.1, 1.1])
                end
                xticks(diff_locs)
            end
        end

        %% Model fitting

        % general setting for all models
        model.num_runs        = num_runs;
        model.x               = (-512:1:512) * deg_per_px;
        model.sA              = sA;
        model.sV              = model.sA;
        model.num_rep         = num_rep;
        model.num_SD          = 5;
        model.numBins_A       = 15;
        model.numBins_V       = 15;
        model.modality        = {'A','V'};
        model.num_rep         = num_rep;
        model.strategy_loc    = 'MA';

        % set OPTIONS to tell bads that my objective function is noisy
        OPTIONS.UncertaintyHandling                            = 1;

        % fit by the corresponding model
        data.org_resp         = org_resp;
        data.org_conf         = org_conf;
        data.sigM             = 1.36; % emperical motor noise averaged from first four participants

        % fit by generating model

        currModel = str2func('nllBimodal');

        % switch confidence strategies
        model.strategy_conf         = ds_conf{d};

        % initiate
        model.mode                  = 'initiate';
        Val = currModel([], model, data);

        % optimize
        model.mode                  = 'optimize';
        NLL                         = NaN(1, model.num_runs);
        estP                        = NaN(model.num_runs, Val.num_para);

%         p = [1 3 0.01025390625 3.336669921875 2.504638671875 20 0.16669921875 33.0859835100118] ;
%         test = currModel(p, model, data);

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
        saveModel{d}.estP(i,:,:) = estP;
        saveModel{d}.NLL(i,:) = NLL;
        saveModel{d}.bestP(i,:) = bestP;
        saveModel{d}.minNLL(i)= minNLL;

        % predict using the best-fitting parameter
        model.mode            = 'predict';
        tmpPred               = currModel(bestP, model, data);
        tmpPred.bestP         = bestP;
        pred{i, d}         = tmpPred;

    end

    save(flnm)

end

%% Plot parameters (predicted vs. ground-truth)
%
% fn                    = fieldnames(GT);
% num_para              = numel(fn);
%
% for d                = 1:num_model
%
%     figure;
%     set(gcf, 'Position', get(0, 'Screensize'));
%     t                     = tiledlayout(2, 4);
%     title(t, sprintf('%s, rep: %i', ds_conf{dd}, num_rep),'FontSize',15);
%
%     % Loop through each parameter
%     for jj                = 1:num_para
%
%         nexttile;
%         hold on
%
%         % Scatter plot of the i-th predicted parameters against its ground-truth
%         scatter(samples(:,jj)', saveModel{d}.bestP(:,jj), 'k','filled');
%         xlabel('Ground Truth');
%         ylabel('Predicted');
%         axis square; % Make the plot square
%         axis equal
%
%         % plot identity line
%         minVal = min([samples(:,jj)', saveModel{d}.bestP(:,jj)]);
%         maxVal = max([samples(:,jj)', saveModel{d}.bestP(:,jj)]);
%         plot([minVal, maxVal], [minVal, maxVal], '--', 'LineWidth', 1);
%
%         % Calculate the Pearson correlation coefficient and p-value
%         [R, P]                = corrcoef(samples(:,jj)', saveModel{d}.bestP(dd,:,jj));
%
%         % Extract the correlation coefficient and p-value
%         r                     = R(1,2);
%         p_value               = P(1,2);
%
%         % Label the r and p in the title
%         title(sprintf('Param %s: r= %.2f, p=%.3f', fn{jj}, r, p_value));
%     end
%
% end
%

