clear; close all;

%% key model recovery parameters

num_rep               = 12;
num_runs              = 20;
num_sample            = 1;
checkFakeData         = true;

%% Manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, 'v1_parameter_recovery');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV1')))

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

    % simulated model info
    ds_loc                = {'Model averaging'};
    ds_conf               = {'Optimal','Suboptimal','Heuristic'};
    num_model             = numel(ds_conf);
    cue_label             = {'Post-cue: A','Post-cue: V'};
    num_cue               = numel(cue_label);
    rel_label             = {'High visual reliability','Low visual reliability'};

    %% Free parameters

    % preallocate
    saveModel             = cell(1, num_model); % for each simulated model
    pred                  = cell(num_sample, num_model, num_model);
    [org_resp, org_conf, org_var]  = deal(NaN(numel(ds_conf),numel(sA), numel(sV), num_cue, numel(rel_label), num_rep));

    for d = 1:num_model

        GT.aA                 = [ 1.1,   1.5]; % scale
        GT.bA                 = [  -1,     1]; % intercept in degree
        GT.sigV1              = [   1,   1.2]; % degree
        GT.delta_sigA         = [ 2, 3]; % degree
        GT.delta_sigV2        = [ 3, 5]; % degree
        GT.sigP               = [  10,    20]; % degree
        GT.pC1                = [ 0.6,   0.7]; % weight
        % calculate criterion based on noise
        [c_lb, ~]= getCriterionRange(GT.sigV1(1), (GT.sigV1(1)+GT.delta_sigA(1)), (GT.sigV1(1)+GT.delta_sigA(1)), GT.sigP(1), ds_conf{d});
        [~, c_ub]= getCriterionRange(GT.sigV1(2), (GT.sigV1(2)+GT.delta_sigA(2)), (GT.sigV1(2)+GT.delta_sigA(2)), GT.sigP(2), ds_conf{d});
        GT.criterion = [c_lb*2, c_ub];

        fn                    = fieldnames(GT);
        num_para              = numel(fn);
        [lb, ub]              = deal(NaN(1,  num_para));
        for k                 = 1:numel(fn)
            lb(:,k)               = GT.(fn{k})(1);
            ub(:,k)               = GT.(fn{k})(2);
        end

        % randomly sample from the range
        samples(d,:,:)               = lb + (ub - lb) .* rand(num_sample, num_para);

        %% Simulate data

        saveModel{d}.estP     = NaN(num_sample, num_runs, num_para);
        saveModel{d}.NLL      = NaN(num_sample, num_runs);
        saveModel{d}.bestP    = NaN(num_sample, num_para);
        saveModel{d}.minNLL   = NaN(num_sample, 1);

        for i                 = 1:num_sample

            disp(i)
            aA                    = samples(d, i, 1);
            bA                    = samples(d, i, 2);
            sigV1                 = samples(d, i, 3);
            sigA                  = samples(d, i, 4) + sigV1;
            sigV2                 = samples(d, i, 5) + sigV1;
            sigVs                 = [sigV1, sigV2];
            sigP                  = samples(d, i, 6);
            pCommon               = samples(d, i, 7);
            criterion             = samples(d, i, 8);

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

            d_resp                = squeeze(loc(:,d,:,:,:));
            org_resp(d,:,:,:,:,:) = reshape(d_resp, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

            d_conf                = squeeze(conf(:,d,:,:,:));
            org_conf(d,:,:,:,:,:) = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

            d_var                = squeeze(variance(:,d,:,:,:));
            org_var(d,:,:,:,:,:) = reshape(d_var, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

        end
    end

    %% check fake data
    if checkFakeData

        for d = 1:2%num_model

            %{diff} cue x reliability x rep
            [conf_by_diff, all_diffs] = org_by_diffs(squeeze(org_conf(d,:,:,:,:,:)), sA);

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
                    ln(rel) = plot(diff_locs, p_conf, 'Color',clt(rel+1,:));
                    patch([diff_locs, fliplr(diff_locs)], ...
                        [p_conf - se_conf, fliplr(p_conf + se_conf)], ...
                        clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1)

                    %                         ylim([-0.1, 1.1])
                end
                xticks(diff_locs)
            end
            legend(ln,rel_label, 'Location', 'eastoutside');

        end
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

        % fit by the same model

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

fn                    = fieldnames(GT);
num_para              = numel(fn);

for dd                = 1:num_model

    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    t                     = tiledlayout(2, 4);
    title(t, sprintf('%s, rep: %i', ds_conf{dd}, num_rep),'FontSize',15);

    % Loop through each parameter
    for jj                = 1:num_para

        nexttile;
        hold on

        % Scatter plot of the i-th predicted parameters against its ground-truth
        scatter(samples(:,jj)', saveModel{d}.bestP(:,jj), 'k','filled');
        xlabel('Ground Truth');
        ylabel('Predicted');
        axis square; % Make the plot square
        axis equal

        % plot identity line
        minVal = min([samples(:,jj)', saveModel{d}.bestP(:,jj)]);
        maxVal = max([samples(:,jj)', saveModel{d}.bestP(:,jj)]);
        plot([minVal, maxVal], [minVal, maxVal], '--', 'LineWidth', 1);

        % Calculate the Pearson correlation coefficient and p-value
        [R, P]                = corrcoef(samples(:,jj)', saveModel{d}.bestP(:,jj));

        % Extract the correlation coefficient and p-value
        r                     = R(1,2);
        p_value               = P(1,2);

        % Label the r and p in the title
        title(sprintf('Param %s: r= %.2f, p=%.3f', fn{jj}, r, p_value));
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

figure; hold on
set(gcf, 'Position',[0,0,270,250]);
% axis equal

imagesc(CM);
colormap('bone')
xticks(1:num_model)
yticks(1:num_model)
xticklabels(ds_conf)
yticklabels(ds_conf)
xlabel('Fit Model');
ylabel('Simulated Model');

[num_rows, num_cols] = size(CM);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM(row, col);
        % Choose text color for better contrast
        textColor = 'w'; % default black
        if val > 0.5
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end


