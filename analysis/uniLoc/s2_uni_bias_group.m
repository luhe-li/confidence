clear; clc; close all;

%% set up

sub_slc = 1;
ses_slc = [1,3];

num_sub = 1;
ses_labels = {'-A','-V1','-V2'};
num_ses = 3;

save_fig = 1;

num_rep = 20;
num_loc = 8;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% clean


[estMu, respSD, sd] = deal(NaN(num_sub, num_ses, 8));

for i = 1:numel(sub_slc)

    sub = sub_slc(i);

    for j = 1:numel(ses_slc)

        ses = ses_slc(j);
        ses_label = ses_labels{ses};
     
        load(sprintf('uniLoc_sub%i_ses%s', sub, ses_label))

        % sort by level
        locRep = reshape([sortedResp(1:end).target_deg],[num_rep,num_loc]);
        loc = locRep(1,:);
        temp_resp = reshape([sortedResp(1:end).response_deg],[num_rep,num_loc]);

        resp(i, j, :, :) = temp_resp; % subject, session(aud, v1, v2), location, rep

        % estimate bias and variance
        respMu(i, j, :) = mean(temp_resp, 1);
        respSD(i, j, :) = std(temp_resp, [], 1);

        % overall variance
        respVar(i, j) = mean(respSD);

    end

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
e = errorbar(loc, squeeze(estMu(sub, s, :)),squeeze(respSD(sub, s, :)),'LineWidth',lw);
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

%% plot variance