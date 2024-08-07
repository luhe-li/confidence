clear; close all; rng('shuffle');

sub_slc = 13;
ses_slc = 1:2; % bimodal sessions

models = {'Heuristic','Suboptimal','Optimal'};

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
result_dir            = fullfile(cur_dir, 'modelFit');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(result_dir))

%% plot set up

lw = 1;
pred_lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% predict data

% load best fitting parameters
flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
files = dir(fullfile(result_dir, flnm));
load(files(end).name);

% winning model?
[~, d] = min([saveConfModel{1}.minNLL, saveConfModel{2}.minNLL, saveConfModel{3}.minNLL]);

% use the best-fitting parameter and winning model
p = saveConfModel{d}.bestP;
p = [p(1), 0.3, 0.2, 0.2, p(5), p(6), p(7)];
model.mode = 'predict';
model.model_slc             = d;
model.strategy_conf         = models{d};

% increase trial number for prediction
model.uni_nrep = 1e3;
model.bi_nrep = 1e3;
pred = nllUniBiConf(p, model, data);
% disp(p)

%% plot bimodal localization response as a function of stimulus location

%TO-DO: use predicted unimodal data

bi_resp = pred.bi_loc;
aud_locs = model.bi_sA;
remapped_vis_locs = model.bi_sV;
raw_diff = unique(aud_locs - aud_locs');
figure; set(gcf, 'Position',[10 10 1200 1000])
plotInd = 1;
for i = 1:2
    for j = 1:2
        subplot(2,2,plotInd)
        cue = i;
        reliability = j;
        plot_spread_VE(bi_resp,aud_locs,raw_diff,remapped_vis_locs,cue,reliability)
        plotInd = plotInd + 1;
    end
end
sgtitle(sprintf('Sub%i, best-fitting model: %s', sub_slc, models{d}), 'FontSize', titleSZ)

%% first panel: plot unimodal localization

figure; 
set(gcf, 'Position',[10 10 1500 800])

subplot(2,3,[1,4])
title('Unimodal localization', 'FontSize',titleSZ)
set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
xlabel('Stimulus location', 'FontSize', titleSZ)
ylabel('Localization','FontSize', titleSZ)
hold on

% sA = sV for uni task
lim = 15;

for j = 1:3

    avg_loc = mean(data.org_uni_loc(j,:,:),3);
    sd_loc = std(data.org_uni_loc(j,:,:), [], 3);

    % data
    plot(model.uni_sA, avg_loc, 'o',...
        'Color', clt(j,:), 'MarkerFaceColor', clt(j,:), 'MarkerSize', 10);
    patch([model.uni_sA, fliplr(model.uni_sA)], ...
        [avg_loc - sd_loc, fliplr(avg_loc + sd_loc)], ...
                clt(j,:),'EdgeColor','none','FaceAlpha',0.1, ...
        'HandleVisibility', 'off');

    % prediction
    plot(model.uni_sA, mean(pred.uni_loc(j,:,:),3),'-','LineWidth',pred_lw,'Color',clt(j,:));

end
plot([-lim, lim], [-lim*1.2, lim*1.2],'k--','LineWidth',lw)

%% second panel: plot unimodal localization

pred.biExpInfo = data.biExpInfo;

[pred_mean_conf, pred_std_mean_conf, pred_uni_pconf, ~,...
    pred_mean_ve, pred_std_ve, ~] = analyze_data(sub_slc, ses_slc, pred);

% analyze real data
[mean_conf, std_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc);

cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High visual reliability','Low visual reliability'};
num_rel               = numel(rel_label);

sgtitle(sprintf('Sub%i, best-fitting model: %s', sub_slc, models{d}),'FontSize',titleSZ)
sp_idx = [2,5];
for cue = 1:num_cue

    subplot(2,3,sp_idx(cue))
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')

    title(cue_label{cue})
    hold on
    xticks(round(raw_diff))
    xtickangle(45)
    ylabel('Shift of localization');

    for rel = 1:num_rel

        % plot data
        plot(raw_diff, squeeze(mean_ve(:, cue, rel)), 'o',...
              'Color', clt(rel+1,:), 'MarkerFaceColor', clt(rel+1,:), 'MarkerSize', 6);
        patch([raw_diff, fliplr(raw_diff)], ...
            [mean_ve(:, cue, rel)' - std_ve(:,cue, rel)', ...
            fliplr(mean_ve(:, cue, rel)' + std_ve(:,cue, rel)')],...
            clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
            'HandleVisibility', 'off');

        % plot prediction
        l(rel) = plot(raw_diff, squeeze(pred_mean_ve(:, cue, rel)), '-','Color',clt(rel+1,:),'LineWidth',pred_lw);

    end

    yline(0,'--')
    if cue == 1
        plot(raw_diff, raw_diff,'k--')
    else
        plot(raw_diff, -raw_diff,'k--')
%         legend([l(:)],rel_label)
        xlabel('Audiovisual discrepancy (V-A, deg)');
    end
end

%% plot confidence report

for cue = 1:num_cue

    subplot(2,3,cue*3)
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')

    title(cue_label{cue})
    hold on
    xticks(round(abs_diff))
    ylim([1 4])
    xlim([min(abs_diff), max(abs_diff)])

    for rel = 1:num_rel

        % plot data
        % Plot the main data points with filled dots
        plot(abs_diff, squeeze(mean_conf(:, cue, rel)), 'o', ...
            'Color', clt(rel+1,:), 'MarkerFaceColor', clt(rel+1,:), 'MarkerSize', 6);

        patch([abs_diff, fliplr(abs_diff)], ...
            [mean_conf(:, cue, rel)' - std_conf(:,cue, rel)', ...
            fliplr(mean_conf(:, cue, rel)' + std_conf(:,cue, rel)')],...
            clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
            'HandleVisibility', 'off');

        % plot prediction
        l(rel) = plot(abs_diff, squeeze(pred_mean_conf(:, cue, rel)), '-', 'Color',clt(rel+1,:),'LineWidth',pred_lw);

        % plot unimodal p_confidence
        if cue == 1
            yline(uni_pconf(1),'--','Color',repmat(0.5,1,3))
        else
            yline(uni_pconf(rel+1),'--','Color',clt(rel+1,:))
            xlabel('Absolute audiovisual discrepancy (deg)');
%             legend([l(:)],rel_label,'Location','best')
        end

    end    
end

saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_conf', sub_slc, d)), 'png')
