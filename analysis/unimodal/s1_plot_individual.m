clear; clc; close all;

%% set up

sub = 'LL';
save_fig = 0;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, mfilename);
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
addpath(genpath(fullfile(project_dir, 'util')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean

[org, exp_info] = org_resp(sub, {'A','V'}, 'uniLoc');
target = org.uni_target;

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

for s = 1:2

    plot(target(s,:), target(s,:), 'k--','LineWidth',lw,'HandleVisibility','off')
    e = errorbar(target(s,:), org.uni_loc_mu(s,:), org.uni_loc_sd(s,:),'Color',clt(s,:),'LineWidth',lw);
    e.CapSize = 0; e.Color = clt(s,:);
    xlim([min(target(s,:))-5, max(target(s,:))+5])
    ylim([min(target(s,:))-5, max(target(s,:))+5])

end

maxlim = 65;
xlim([-maxlim, maxlim])
ylim([-maxlim, maxlim])
title(sprintf('S%i', sub))

xlabel('Stimulus location (cm)')
ylabel('Estimate location (cm)')
legend('A','V')
