clear; close all; rng('shuffle');

sub_slc = 13;
ses_slc = 1:2; % bimodal sessions

models = {'Heuristic','Suboptimal','Optimal'};

str = 'nomuP_pinSigA';

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

%% predict data

% load best fitting parameters
flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
files = dir(fullfile(result_dir, flnm));
load(files(end).name);

% winning model?
[~, d] = min([saveModel{1}.minNLL, saveModel{2}.minNLL, saveModel{3}.minNLL]);

% use the best-fitting parameter and winning model
p = saveModel{d}.bestP;
model.uni_nrep = 1e3;
model.bi_nrep = 1e3;
model.mode = 'predict';
model.model_slc             = d;
model.strategy_conf         = models{d};
pred = nllUniBiLocConf(p, model, data);
disp(p)

%% unimodal

figure; 
set(gcf, 'Position',[10 10, 500 400])

% subplot(2,3,[1,4])
title('Unimodal localization', 'FontSize',titleSZ)
set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
xlabel('Stimulus location', 'FontSize', titleSZ)
ylabel('Localization','FontSize', titleSZ)
hold on

% sA = sV for uni task
lim = 15;
    [data.org_uni_loc, data.org_uni_conf, ~, data.uniExpInfo, ~, ~, data.uni_loc, data.uni_conf, data.coefsA] = org_data(sub_slc,[],'uniLoc');

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
    plot(model.uni_sA, mean(pred.uni_loc(j,:,:),3),'-','LineWidth',2,'Color',clt(j,:));

end
plot([-lim, lim], [-lim*1.2, lim*1.2],'k--','LineWidth',lw)

saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_uniloc_%s', sub_slc, d, str)), 'png')

%% plot VE

% analyze data prediction 
pred.biExpInfo = data.biExpInfo;

[pred_mean_ve, pred_std_ve] = analyze_loc(sub_slc, ses_slc, pred);

% analyze real data
[~, ~, ~, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc);

cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High visual reliability','Low visual reliability'};
num_rel               = numel(rel_label);

figure; hold on
t = tiledlayout(1, 2);
title(t,sprintf('Sub%i, best-fitting model: %s', sub_slc, models{d}))
xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
ylabel(t, 'Shift of localization');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:num_cue

    nexttile
    title(cue_label{cue})
    axis equal
    hold on
    xticks(round(raw_diff))

    for rel = 1:num_rel

        % plot data
        plot(raw_diff, squeeze(mean_ve(:, cue, rel)), 'o',...
              'Color', clt(rel+1,:), 'MarkerFaceColor', clt(rel+1,:), 'MarkerSize', 6);
        patch([raw_diff, fliplr(raw_diff)], ...
            [mean_ve(:, cue, rel)' - std_ve(:,cue, rel)', ...
            fliplr(mean_ve(:, cue, rel)' + std_ve(:,cue, rel)')],...
            clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
            'HandleVisibility', 'off');

%         % Add error bars
%         errorbar(raw_diff, squeeze(mean_ve(:, cue, rel)), squeeze(std_ve(:, cue, rel)), 'o', ...
%             'Color', clt(rel+1,:), 'LineStyle', 'none', 'LineWidth', 2, 'CapSize', 0);

        % plot prediction
        l(rel) = plot(raw_diff, squeeze(pred_mean_ve(:, cue, rel)), '-','Color',clt(rel+1,:),'LineWidth',lw);

    end

    yline(0,'--')
    if cue == 1
        plot(raw_diff, raw_diff,'k--')
    else
        plot(raw_diff, -raw_diff,'k--')
        legend([l(:)],rel_label)
    end
end

saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_loc_%s', sub_slc, d, str)), 'png')

%% plot bimodal localization response as a function of stimulus location

bi_resp = pred.bi_loc;
aud_locs = model.bi_sA;
remapped_vis_locs = model.bi_sV;
raw_diff = unique(aud_locs - aud_locs');
plotInd = 1;
figure; 
set(gcf, 'Position',[10 10, 1500 1000])
for i = 1:2
    for j = 1:2
        subplot(2,2,plotInd)
        cue = i;
        reliability = j;
        plot_spread_VE(bi_resp,aud_locs,raw_diff,remapped_vis_locs,cue,reliability)
        plotInd = plotInd + 1;
    end
end

saveas(gca, fullfile(out_dir, sprintf('sub%i_bestfittingM%i_spread_%s', sub_slc, d, str)), 'png')