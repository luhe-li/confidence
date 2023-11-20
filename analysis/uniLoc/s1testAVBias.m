clear; clc; close all;

%% set up

sub = 1;
ses = '-A';

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean
load(sprintf('uniLoc_sub%i_ses%s', sub, ses))

% basic info
nRep = ExpInfo.nRep;
nLevel = ExpInfo.nLevel;

% sort by level
locRep = reshape([sortedResp(1:end).loc_cm],[ExpInfo.nRep,ExpInfo.nLevel]);
loc = locRep(1,:);
est = reshape([sortedResp(1:end).Response_deg],[ExpInfo.nRep,ExpInfo.nLevel]);

estMu = mean(est, 1);
sdMu = std(est, [], 1);

%% plot

lw = 2;
fontSZ = 15;
titleSZ = 20;
clt = [202,0,32; 5,113,176]./255; % red and blue

figure;hold on
set(gca, 'LineWidth', 2, 'FontSize', 15,'TickDir', 'out')
set(gcf, 'Position',[10 10 500 400])
axis equal

plot(loc, loc,'k--','LineWidth',lw)
e = errorbar(loc, estMu, sdMu,'LineWidth',lw);
e.CapSize = 0; e.Color = clt(2,:);
xlim([min(loc)-5, max(loc)+5])
ylim([min(loc)-5, max(loc)+5])

switch ses

    case '-A'
        xlabel('Auditory stimulus location (cm)')
        ylabel('Visual estimate location (cm)')

end

if save_fig
    flnm = sprintf('sub%i_loc%s',sub, ses);
    saveas(gca, fullfile(out_dir, flnm),'png')
end