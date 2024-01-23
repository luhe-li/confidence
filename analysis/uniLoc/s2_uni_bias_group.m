clear; clc; close all;

%% set up

sub_slc = 4;

% session
ses_labels = {'-A','-V'};

% condition
cond_label = {'A','V-high reliability','V-low reliability'};
num_cond = numel(cond_label);

% fixed parameters
num_loc = 8;
num_rep = 20;

save_fig = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's2Fig');
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% organize data

% load all sessions

for i = 1:numel(sub_slc)

    sub = sub_slc(i);

    for s = 1:2

        ses_label = ses_labels{s};
        load(sprintf('uniLoc_sub%i_ses%s', sub, ses_label))

        if strcmp(ses_label,'-A')

            j = 1; % condition A
            orgResp{j} = sortedResp;

        elseif strcmp(ses_label,'-V')

            j = 2; % condition V1
            orgResp{j} = sortedReli1Resp;

            j = 3; % condition V2
            orgResp{j} = sortedReli2Resp;

        end
    end
end

resp = NaN(numel(sub_slc), num_cond, num_loc, num_rep);
[stim, respMu, respSD] = deal(NaN(numel(sub_slc), num_cond, num_loc));
respVar = NaN(numel(sub_slc), num_cond);

for i = 1:numel(sub_slc)

    sub = sub_slc(i);

    for j = 1:num_cond

        % reshape by stimulus level
        locRep = reshape([orgResp{j}(1:end).target_deg],[num_rep,num_loc]);
        stim(i, j, :) = locRep(1,:);
        temp_resp = reshape([orgResp{j}(1:end).response_deg],[num_rep,num_loc]);

        resp(i, j, :, :) = temp_resp'; % subject, session(aud, v1, v2), location, rep

        % estimate bias and variance
        respMu(i, j, :) = mean(temp_resp, 1);
        respSD(i, j, :) = std(temp_resp, [], 1);

        % overall variance
        respVar(i, j) = mean(respSD(i, j, :));

    end

end

%% plot set up

% figure set up
lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;

clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    251, 154, 153]./255; % light red

%% 1. plot raw responses

stim_level = squeeze(stim(1,1,:));

for i = 1:numel(sub_slc)

    sub = sub_slc(i);

    figure;
    set(gcf, 'Position',[0 0 1400 400])
    t = tiledlayout(1,3);
    xlabel(t,'Stimulus location', 'FontSize', titleSZ)
    ylabel(t,'Estimate','FontSize', titleSZ)

    lim = max(squeeze(resp(i, :, :, :)),[],'all');

    % for each condition
    for j = 1:num_cond

        nexttile(j)
        set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
        axis equal
        hold on
       
        title(cond_label{j})       
        plot([-lim, lim], [-lim, lim],'k--','LineWidth',lw)

        % plot responses for each stimulus location
        for k = 1:numel(stim_level)

            scatter(repmat(stim_level(k), 1, num_rep), squeeze(resp(i, j, k, :)),...
                dotSZ,'MarkerFaceColor',clt(j,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)

        end

    end

end

if save_fig
    flnm = sprintf('sub%i_raw',sub);
    saveas(gca, fullfile(out_dir, flnm),'png')
end

%% 2. plot response mean to check audiovisual bias

for i = 1:numel(sub_slc)

    sub = sub_slc(i);

    figure;
    set(gcf, 'Position',[0 0 500 400])
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
    xlabel('Stimulus location', 'FontSize', titleSZ)
    ylabel('Estimate','FontSize', titleSZ)
    hold on

    lim = max(stim_level,[],'all');

    % for each condition
    for j = 1:num_cond
        
%         plot([-lim, lim], [-lim, lim],'k--','LineWidth',lw)
        e = errorbar(stim_level, squeeze(respMu(i, j, :)),squeeze(respSD(i, j, :)),'LineWidth',lw,'Color',clt(j,:));
        e.CapSize = 0; %e.Color = clt(2,:);

    end

    legend(cond_label,'Location','northwest');
    legend boxoff

end

if save_fig
    flnm = sprintf('sub%i_bias',sub);
    saveas(gca, fullfile(out_dir, flnm),'png')
end

%% 3. plot response variance to check reliability manipulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Make a bar plot of localization response variance for each condition.
% Variance should be averaged across locations. Use the same color code for
% each condition (use `clt`).
% 2. Make a group-average bar plot of localization response variance for
% each condition. Add error bars for across-subject variation.

% YOUR CODES GO HERE %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if save_fig
    flnm = sprintf('sub%i_var',sub);
    saveas(gca, fullfile(out_dir, flnm),'png')
end