clear; clc; close all;

%% set up

sub = 'LL';
save_fig = 0;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, mfilename);
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean

all_ses = {'A','V'};

for s = 1:numel(all_ses)

    load(sprintf('uniLoc_sub-%s_ses-%s', sub, all_ses{s}));

    % basic info
    nRep = ExpInfo.nRep;
    nLevel = ExpInfo.nLevel;

    % modality x location level x repetiton
    loc_rep =  reshape([sortedResp(1:end).target_cm],[nRep, nLevel]); 
    loc(s,:) = loc_rep(1,:); 
    est(s,:,:) = reshape([sortedResp(1:end).response_cm],[nRep, nLevel]);

    % estimate bias and variance
    estMu(s,:) = mean(squeeze(est(s,:,:)), 1);
    sdMu(s,:) = std(squeeze(est(s,:,:)), [], 1);

end

% 
% % sort by level
% locRep = reshape([sortedResp(1:end).target_deg],[nRep,nLevel]);
% loc = locRep(1,:);
% est = reshape([sortedResp(1:end).response_deg],[nRep,nLevel]);
% 
% % estimate bias and variance
% estMu = mean(est, 1);
% sdMu = std(est, [], 1);
% 
% % overall variance
% sd = std(est - locRep, [],"all");

%% plot

% figure set up
lw = 2;
fontSZ = 15;
titleSZ = 20;
clt = [5,113,176; 202,0,32]./255; % red and blue

figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 15,'TickDir', 'out')
set(gcf, 'Position',[10 10 500 400])
hold on
title(sprintf('%s', sub))
axis equal

for s = 1:numel(all_ses)

    plot(loc(s,:), loc(s,:), 'k--','LineWidth',lw)
    e = errorbar(loc(s,:), estMu(s,:), sdMu(s,:),'Color',clt(s,:),'LineWidth',lw);
    e.CapSize = 0; e.Color = clt(s,:);
    xlim([min(loc(s,:))-5, max(loc(s,:))+5])
    ylim([min(loc(s,:))-5, max(loc(s,:))+5])

end

maxlim = 50;
xlim([-maxlim, maxlim])
ylim([-maxlim, maxlim])
title(sprintf('S%i', sub))

switch ses

    case '-A'
        xlabel('Auditory stimulus location (dva)')
        ylabel('Estimate location (dva)')

end

if save_fig
    flnm = sprintf('sub%i_cond%s',sub, ses);
    saveas(gca, fullfile(out_dir, flnm),'png')
end