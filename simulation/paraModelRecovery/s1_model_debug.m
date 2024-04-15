% simulate one set of parameter
% fit by v2: fitAllModels (know inside)

clear; close all;

%% key model recovery parameters

num_rep               = 100;
num_runs              = 4;
num_sample            = 1;
checkFakeData         = 1;

%% Manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
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
    x                     = -screenXdeg /2 : deg_per_px : screenXdeg/2; % the screen deg space

    sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
    sV                    = sA;
    sAV                   = combvec(sA, sV);
    num_s                 = size(sAV,2); % number of conditions

    aud_locs              = sA;
    vis_locs              = sV;
    diffs                 = zeros(length(aud_locs), length(vis_locs));
    for ii                 = 1:length(aud_locs)
        for j                 = 1:length(vis_locs)
            diffs(ii, j)           = aud_locs(ii) - vis_locs(j);
        end
    end
    abs_diff             = unique(abs(diffs))';

    %% Free parameters

    % choose a reasonble set of parameter set. See variable name below.
    %    aA, bA, sigV1, dsigA, dsigV2, sigP,  pCC, sigC, cA, cV
    GT = {[1,  0.1,  1,   1.2,    1.5,   8, 0.57,  0.3, 0.5, 0.5],...% Heuristic
        [1,  0.1,  0.5,   1.5,    1.8,   8,  0.8, 0.5, 0.48, 0.6],...% Suboptimal
        [1,  0.1,  1,   1.5,    1.8,   8,  0.57,  0.5, 0.48, 0.6]}; % Optimal

    % simulated model info
    ds_loc                = {'Model averaging'};
    ds_conf               = {'Heuristic','Suboptimal','Optimal'};
    num_model             = numel(ds_conf);
    cue_label             = {'Post-cue: A','Post-cue: V'};
    num_cue               = numel(cue_label);
    rel_label             = {'High reliability','Low reliability'};

    %% Simulate data

    [fake_data, saveModel, pred] = deal(cell(num_model, num_sample));

    for i = 1:num_sample

        for d = 3%3:-1:1


            % jitter each parameters a little
            j_gt = addRandomJitter(GT{d});

            % assign simulation parameters
            num_para              = length(j_gt);
            aA                    = j_gt(1);
            bA                    = j_gt(2);
            sigV1                 = j_gt(3);
            sigA                  = j_gt(4);
            sigV2                 = j_gt(5);
            sigVs                 = [sigV1, sigV2];
            sigP                  = j_gt(6);
            pCommon               = j_gt(7);
            sigM                  = j_gt(8);
            cA                    = j_gt(9);
            cV                    = j_gt(10);
            lapse                 = 0.02;
            muP                   = 0;

            [loc, conf] = deal(NaN(num_s, num_cue, numel(sigVs), num_rep));
            % num_s: stimulus location combination
            % 3 confidence decision strategies
            % 2 modalities(1 = aud, 2 = vis)
            % 2 visual reliabilities
            % num_rep

            for j                 = 1:num_s
                for v                 = 1:numel(sigVs)
                    [loc(j,:,v,:), conf(j,:,v,:)] = sim_loc_conf(num_rep, sAV(1,j), sAV(2,j),...
                        aA, bA, sigA, sigVs(v), muP, sigP, pCommon, sigM, cA, cV, lapse, d);
                end
            end

            %% organize data

            org_loc               = reshape(loc, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);
            org_conf              = reshape(conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

            data{d,i}.org_loc = org_loc;
            data{d,i}.org_conf = org_conf;
            data{d,i}.gt = j_gt;
            data{d,i}.sigMotor         = 1.36; % emperical motor noise averaged from first four participants

            %% check fake data
            if checkFakeData

                %% check localization data

                uni_loc = zeros(size(org_loc));

                loc_a = repmat(sA',[1,numel(sV)]);
                loc_v = repmat(sV,[numel(sA),1]);

                uni_loc(:,:,1,1,:) = repmat(loc_a, [1, 1, 1, 1, num_rep]);
                uni_loc(:,:,2,1,:) = repmat(loc_v, [1, 1, 1, 1, num_rep]);

                uni_loc(:,:,1,2,:) = uni_loc(:,:,1,1,:);
                uni_loc(:,:,2,2,:) = uni_loc(:,:,2,1,:);

                % loc at uni minus loc at bi
                ve =  mean(org_loc,5) - mean(uni_loc, 5);

                % diff x cue x reliability
                [ve_by_raw_diff, raw_diff] = org_by_raw_diffs_4D(ve, sA);

                % assume participants localized perfectly in the unisensory
                % condition
                figure; hold on
                t = tiledlayout(2, 1);
                title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
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
                    xticks(raw_diff)
                    yline(0,'--')
                    if cue == 1
                        plot(raw_diff, raw_diff,'k--')
                    else
                        plot(raw_diff, -raw_diff,'k--')
                    end
                end

                %% check confidence data

                % organize localization data: {diff} cue x reliability x rep
                [conf_by_diff, all_diffs] = org_by_diffs(org_conf, sA);

                figure; hold on
                t = tiledlayout(2, 1);
                title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
                xlabel(t, 'Absolute audiovisual discrepancy (deg)');
                ylabel(t, 'Proportion of confidence report');
                t.TileSpacing = 'compact';
                t.Padding = 'compact';

                for cue = 1:num_cue
                    nexttile
                    title(cue_label{cue})
                    hold on
                    for rel = 1: numel(sigVs)
                        [p_conf, se_conf] = deal(NaN(1, numel(abs_diff)));
                        for diff = 1:numel(abs_diff)
                            i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
                            p = sum(i_conf)/numel(i_conf);
                            p_conf(diff) = p;
                            se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
                        end
                        plot(abs_diff, p_conf, 'Color',clt(rel+1,:));
                        ylim([-0.1, 1.1])
                    end
                    xticks(abs_diff)
                end
            end

            %% Model fitting

            % general setting for all models
            model.num_runs        = num_runs;
            model.num_sec         = 10; % number of samples in the parameter space, must be larger than num_runs
            model.x               = (-512:1:512) * deg_per_px;
            model.sA              = sA;
            model.sV              = model.sA;
            model.num_rep         = num_rep;
            model.num_SD          = 5;
            model.numBins_A       = 15;
            model.numBins_V       = 15;
            model.modality        = {'A','V'};
            model.strategy_loc    = 'MA';

            % set OPTIONS to tell bads that my objective function is noisy
            %         OPTIONS.UncertaintyHandling     = 0;
            OPTIONS.TolMesh = 1e-3;

            % fit by the corresponding model
            data.org_resp         = org_loc;
            data.org_conf         = org_conf;
            data.sigMotor         = 1.36; % emperical motor noise averaged from first four participants

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

            % % test using ground truth parameters
            %             p = GT{d};
            %             test = currModel(p, model, data);

            parfor n              = 1:model.num_runs

                tempModel             = model;
                tempVal               = Val;
                tempFunc              = currModel;

                [estP(n,:),NLL(n),~,~,traj(n)]    = bads(@(p) tempFunc(p, model, data),...
                    Val.init(n,:), Val.lb,...
                    Val.ub, Val.plb, Val.pub, [], OPTIONS);

                disp(estP(n,:))

            end

            % find the parameter with the least NLL
            [minNLL, best_idx]    = min(NLL);
            bestP                 = estP(best_idx, :);

            % save all fitting results
            saveModel{d,i}.estP = estP;
            saveModel{d,i}.NLL = NLL;
            saveModel{d,i}.bestP = bestP;
            saveModel{d,i}.minNLL = minNLL;

            % predict using the best-fitting parameter
            model.mode            = 'predict';
            tmpPred               = currModel(bestP, model, data);
            pred{d,i}             = tmpPred;

        end
    end

    save(flnm, 'saveModel','pred','data')

end

%% Plot parameters (predicted vs. ground-truth)

% load data

% set other info
model.num_runs        = 10;
model.num_sec         = 10; 
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
num_rep               = 100;
num_sample            = 100;

currModel = str2func('nllBimodal');
model.mode                  = 'initiate';
init = currModel([], model, data);

num_para              = init.num_para;
paraID                = init.paraID;

for d                = 1:num_model

    figure;
    set(gcf, 'Position', [1 188 1920 789]);
    t  = tiledlayout(2, 5);
    title(t, sprintf('%s, rep: %i', ds_conf{d}, num_rep),'FontSize',15);
    xlabel(t, 'Ground Truth');
    ylabel(t, 'Predicted');

    % Loop through each parameter
    for jj                = 1:num_para

        nexttile;
        hold on

        for i = 1:num_sample

            % extract
            gt = data{d,i}.gt(jj);
            pr = saveModel{d,i}.bestP(jj);

            % save in a matrix for correlation
            all_gt(i, jj) = gt;
            all_pr(i, jj) = pr;

            if i == num_sample

                % Scatter plot of the i-th predicted parameters against its ground-truth
                scatter(all_gt(:, jj), all_pr(:, jj), 50,'k','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                axis square; % Make the plot square
                axis equal

                % identity line
                plot([init.lb(jj), init.ub(jj)], [init.lb(jj), init.ub(jj)], 'k--', 'LineWidth', 1);

                % Calculate the Pearson correlation coefficient and p-value
                [R, P]                = corrcoef(all_gt(:, jj), all_pr(:, jj));

                % Extract the correlation coefficient and p-value
                r                     = R(1,2);
                p_value               = P(1,2);

                % Label the r and p in the title
                xlabel(sprintf('%s: r= %.2f, p=%.3f', paraID{jj}, r, p_value));

            end

        end

    end

    saveas(gca, fullfile(pwd, 's1_model_debug', sprintf('recovery_model-%s_rep-%i', ds_conf{d}, num_rep)), 'png')
    
end


%% function section

function newValue = addRandomJitter(originalValue)
% Calculate 10% of the original value
jitterPercent = 0.01; % 10%
jitterAmount = originalValue * jitterPercent;
newValue = NaN(size(originalValue));
% Randomly choose to add or subtract the jitter
randInd = randi([0 1], size(originalValue));
% Add jitter
newValue(~~randInd) = originalValue(~~randInd) + jitterAmount(~~randInd);
% Subtract jitter
newValue(~randInd) = originalValue(~randInd) - jitterAmount(~randInd);

end
