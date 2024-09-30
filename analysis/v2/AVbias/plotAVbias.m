% Compare whether audiovisual bias changes each session (1 unimodal, 2-3
% bimodal)

clear; clc; close all;

sub_slc = 1;
ses_slc = 1:2; % num of bimodal session
nses = numel(ses_slc);

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
out_dir      = fullfile(pwd, mfilename);
addpath(genpath(fullfile(project_dir,'func')))
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load and organize data

% condition (A,V1,V2) x loc (4) x rep
uni_resp = org_data(sub_slc,[],'uniLoc');

% location x response modality
uni_mu(:, 1) = squeeze(mean(uni_resp(1,:,:), 3));
uni_mu(:, 2) = squeeze(mean(mean(uni_resp(2:3,:,:), 1),3));

uni_sd(:, 1) = std(squeeze(uni_resp(1,:,:)),[],2);
uni_sd(:, 2) = std(squeeze(mean(uni_resp(2:3,:,:), 1)),[],2);

for ses = 1:nses
    bi_resp = org_data(sub_slc,ses,'biLoc');

    % only look at congruent trials
    for loc = 1:4

        % modality index: 1 = a, 2 = v
        for mm = 1:2

            % session x location x response modality
            curr_resp = squeeze(bi_resp(loc, loc, mm, :, :));
            bi_mu(ses, loc, mm) = mean(curr_resp,'all');
            bi_sd(ses, loc, mm)  = std(curr_resp,[],'all');

        end
    end
end

%% plot

lw = 1;
fontSZ = 15;

figure; hold on
set(gcf, 'Position', [0,0,1100,300]); 

% unimodal
subplot(1,1+nses, 1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
axis equal; axis square;

errorbar(uni_mu(:,2), uni_mu(:,1), uni_sd(:,1),'k','LineWidth',lw,'CapSize',0);
errorbar(uni_mu(:,2), uni_mu(:,1), uni_sd(:,2),'k','horizontal','LineWidth',lw,'CapSize',0);
xlabel('Visual localization')
ylabel('Auditory localization')
title('Unimodal')

% Determine the limits for the identity line
minLimit = -30;
maxLimit = 30;

% Plot the identity line
plot([minLimit, maxLimit], [minLimit, maxLimit], 'k--');
xlim([minLimit, maxLimit]);
ylim([minLimit, maxLimit]);
xticks([-20, 0, 20])
yticks([-20, 0, 20])

for ses = 1:nses

    subplot(1,1+nses, 1+ses); hold on
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
    axis equal; axis square;

    ses_bi_mu = squeeze(bi_mu(ses,:,:));
    ses_bi_sd = squeeze(bi_sd(ses,:,:));

    errorbar(ses_bi_mu(:,2), ses_bi_mu(:,1), ses_bi_sd(:,1),'k','LineWidth',lw,'CapSize',0);
    errorbar(ses_bi_mu(:,2), ses_bi_mu(:,1), ses_bi_sd(:,2),'k','horizontal','LineWidth',lw,'CapSize',0);
    xlabel('Visual localization')
    ylabel('Auditory localization')
    title(sprintf('Bimodal ses-%i',ses))

    % Plot the identity line
    plot([minLimit, maxLimit], [minLimit, maxLimit], 'k--');
    xlim([minLimit, maxLimit]);
    ylim([minLimit, maxLimit]);

    xticks([-20, 0, 20])
    yticks([-20, 0, 20])

end

sgtitle(sprintf('sub%i, bimodal corrected', sub_slc))

flnm = sprintf('avbias_check_sub%i_ses%i-%i', sub_slc, min(ses_slc), max(ses_slc));
saveas(gca, fullfile(out_dir, flnm),'png')

