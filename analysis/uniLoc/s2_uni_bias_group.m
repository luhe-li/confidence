clear; clc; close all;

%% set up

sub = 4;
sess = {'-A','-V1','-V2'};
save_fig = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean

[estMu, sdMu, sd] = deal(NaN(15, 3, 8));

for s = 1%[1,3]

    ses = sess{s};
    load(sprintf('uniLoc_sub%i_ses%s', sub, ses))

    % basic info
    nRep = ExpInfo.nRep;
    nLevel = ExpInfo.nLevel;

    % sort by level
    locRep = reshape([sortedResp(1:end).target_deg],[nRep,nLevel]);
    loc = locRep(1,:);
    temp_est = reshape([sortedResp(1:end).response_deg],[nRep,nLevel]);

    est(sub, s, :, :) = temp_est; % subject, session(aud, v1, v2), location, rep

    % estimate bias and variance
    estMu(sub, s, :) = mean(temp_est, 1);
    sdMu(sub, s, :) = std(temp_est, [], 1);

    % overall variance
%     sd(sub, s, :)= std(temp_est - locRep, [],"all");

end



%% plot bias

% figure set up
lw = 2;
fontSZ = 15;
titleSZ = 20;
clt = [202,0,32; 5,113,176]./255; % red and blue

figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 15,'TickDir', 'out')
set(gcf, 'Position',[10 10 500 400])
axis equal

plot(loc, loc,'k--','LineWidth',lw)
e = errorbar(loc, squeeze(estMu(sub, s, :)),squeeze(sdMu(sub, s, :)),'LineWidth',lw);
e.CapSize = 0; e.Color = clt(2,:);
xlim([min(loc)-5, max(loc)+5])
ylim([min(loc)-5, max(loc)+5])
title(sprintf('S%i', sub))

switch ses

    case '-A'
        xlabel('Auditory stimulus location (dva)')
        ylabel('Visual estimate location (dva)')

end

if save_fig
    flnm = sprintf('sub%i_loc%s',sub, ses);
    saveas(gca, fullfile(out_dir, flnm),'png')
end

%% plot variance