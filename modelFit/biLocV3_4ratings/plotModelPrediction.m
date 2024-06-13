clear; close all; rng('shuffle');

sub_slc = 13;
ses_slc = 1:3; % bimodal sessions

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

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;
%     251, 154, 153]./255; % light red

%% predict data

% load best fitting parameters
flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
files = dir(fullfile(result_dir, flnm));
load(files(end).name);

% winning model?
[~, d] = min([saveConfModel{1}.minNLL, saveConfModel{2}.minNLL, saveConfModel{3}.minNLL]);
d=2;

% use the best-fitting parameter and winning model
p = saveConfModel{d}.bestP;
model.mode = 'predict';
model.model_slc             = d;
model.strategy_conf         = models{d};
pred = nllLocConf(p, model, data);
disp(p)

% %% plot bimodal localization response as a function of stimulus location
% 
% bi_resp = pred.bi_loc;
% aud_locs = model.bi_sA;
% remapped_vis_locs = model.bi_sV;
% raw_diff = unique(aud_locs - aud_locs');
% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% plotInd = 1;
% for i = 1:2
%     for j = 1:2
%         subplot(2,2,plotInd)
%         cue = i;
%         reliability = j;
%         plot_spread_VE(bi_resp,aud_locs,raw_diff,remapped_vis_locs,cue,reliability)
%         plotInd = plotInd + 1;
%     end
% end
%% plot VE

% analyze prediction 
pred.biExpInfo = data.biExpInfo;
[pred_mean_conf, pred_diff] = analyze_pred(pred);

% [pred_mean_conf, pred_std_mean_conf, pred_uni_pconf, ~,...
%     pred_mean_ve, pred_std_ve, ~] = analyze_data(sub_slc, ses_slc, pred, false);

% analyze real data
[mean_conf, std_conf, uni_pconf, data_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc);

cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High visual reliability','Low visual reliability'};
num_rel               = numel(rel_label);
% 
% figure; hold on
% t = tiledlayout(1, 2);
% title(t,sprintf('Sub%i, best-fitting model: %s', sub_slc, models{d}))
% xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
% ylabel(t, 'Shift of localization');
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% 
% for cue = 1:num_cue
% 
%     nexttile
%     title(cue_label{cue})
%     axis equal
%     hold on
%     xticks(round(raw_diff))
% 
%     for rel = 1:num_rel
% 
%         % plot data
%         plot(raw_diff, squeeze(mean_ve(:, cue, rel)), 'o',...
%               'Color', clt(rel+1,:), 'MarkerFaceColor', clt(rel+1,:), 'MarkerSize', 6);
%         patch([raw_diff, fliplr(raw_diff)], ...
%             [mean_ve(:, cue, rel)' - std_ve(:,cue, rel)', ...
%             fliplr(mean_ve(:, cue, rel)' + std_ve(:,cue, rel)')],...
%             clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
%             'HandleVisibility', 'off');
% 
% %         % Add error bars
% %         errorbar(raw_diff, squeeze(mean_ve(:, cue, rel)), squeeze(std_ve(:, cue, rel)), 'o', ...
% %             'Color', clt(rel+1,:), 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 0);
% 
%         % plot prediction
%         l(rel) = plot(raw_diff, squeeze(pred_mean_ve(:, cue, rel)), '-','Color',clt(rel+1,:),'LineWidth',lw);
% 
%     end
% 
%     yline(0,'--')
%     if cue == 1
%         plot(raw_diff, raw_diff,'k--')
%     else
%         plot(raw_diff, -raw_diff,'k--')
%         legend([l(:)],rel_label)
%     end
% end
% 
% saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_loc', sub_slc, d)), 'png')

%% plot confidence report

figure; hold on
t = tiledlayout(2, 1);
title(t,sprintf('Sub%i, best-fitting model: %s', sub_slc, models{d}))
xlabel(t, 'Absolute audiovisual discrepancy (deg)');
ylabel(t, 'Proportion of confidence report');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:num_cue

    nexttile
    title(cue_label{cue})
    hold on
    xticks(round(data_diff))
    ylim([1 4])
    xlim([min(data_diff), max(data_diff)])

    for rel = 1:num_rel

        % plot data
        % Plot the main data points with filled dots
        plot(data_diff, squeeze(mean_conf(:, cue, rel)), 'o', ...
            'Color', clt(rel+1,:), 'MarkerFaceColor', clt(rel+1,:), 'MarkerSize', 6);

        patch([data_diff, fliplr(data_diff)], ...
            [mean_conf(:, cue, rel)' - std_conf(:,cue, rel)', ...
            fliplr(mean_conf(:, cue, rel)' + std_conf(:,cue, rel)')],...
            clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
            'HandleVisibility', 'off');

        % plot prediction
        l(rel) = plot(pred_diff, squeeze(pred_mean_conf(:, cue, rel)), '-', 'Color',clt(rel+1,:),'LineWidth',lw);

        % plot unimodal p_confidence
        if cue == 1
            yline(uni_pconf(1),'--','Color',repmat(0.5,1,3))
        else
            yline(uni_pconf(rel+1),'--','Color',clt(rel+1,:))
            legend([l(:)],rel_label,'Location','best')
        end

    end    
end

saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_conf', sub_slc, d)), 'png')