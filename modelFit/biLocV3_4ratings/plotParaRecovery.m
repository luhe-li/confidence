clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
results_dir           = fullfile(cur_dir, 'paraRecovery');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% load fitting results

num_rep = 46; % experimental repetition per condition
num_run = 7; % run of fits
num_sample = 100; % sample of GT
flnm = sprintf('recoveryResults_rep%i_sample%i_run%i.mat', num_rep, num_sample, num_run);

if ~exist(fullfile(results_dir, flnm),'file') 
    sprintf('Parameter recovery results don''t exist')
else
    load(fullfile(results_dir, flnm));
end

%% plot ground-truth against loc stage & loc+conf stage for each model

ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
paraID = saveConfModel{1,1}.paraInfo.paraID;
num_para = numel(paraID);

% take lb and ub from the loc+conf model
lb = saveConfModel{1,1}.paraInfo.lb;
ub = saveConfModel{1,1}.paraInfo.ub;

for d                = 1:num_model

    figure;
    set(gcf, 'Position', [1 188 1500 789]);
    t  = tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('%s, rep: %i', ds_conf{d}, num_rep),'FontSize',20);
%     xlabel(t, 'Ground Truth','FontSize',15);
%     ylabel(t, 'Predicted','FontSize',15);

    % Loop through each parameter
    for jj                = 1:num_para

        nexttile;
        hold on

        for i = 1:num_sample

            % extract
            if jj <= 7
                pr1 = saveLocModel{d,i}.bestP(jj);
                all_pr1(i, jj) = pr1;
            end
            gt = saveData{d,i}.gt(jj);
            pr2 = saveConfModel{d,i}.bestP(jj);

            % save in a matrix for correlation
            all_gt(i, jj) = gt;
            all_pr2(i, jj) = pr2;

            if i == num_sample

                if jj <= 7
                    % Scatter plot of the i-th predicted parameters from loc-model against its ground-truth
                    s1 = scatter(all_gt(:, jj), all_pr1(:, jj), 50,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); 
                end

                % loc+conf model
                s2 = scatter(all_gt(:, jj), all_pr2(:, jj), 50,'MarkerEdgeColor','k','MarkerFaceColor','none');
                axis square;
                axis equal

                % identity line
                ax = gca;
                x_limits = ax.XLim;
                y_limits = ax.YLim;
                line_min = min([x_limits y_limits]);
                line_max = max([x_limits y_limits]);
                plot([line_min line_max], [line_min line_max], 'k--', 'LineWidth', 1);

                % Calculate the Pearson correlation coefficient and p-value
                [R, P]                = corrcoef(all_gt(:, jj), all_pr2(:, jj));

                % Extract the correlation coefficient and p-value
                r                     = R(1,2);
                p_value               = P(1,2);

                % Label the r and p in the title
                title(sprintf('%s: r= %.2f, p=%.3f', paraID{jj}, r, p_value));
                xlabel('Ground-truth')
                ylabel('Prediction')

            end

        end

    end

    lgd = legend([s1, s2],{'Loc only', 'Loc & Conf'});
    lgd.Layout.Tile = 'east';

    saveas(gca, fullfile(out_dir, sprintf('recovery_model-%s_rep-%i', ds_conf{d}, num_rep)), 'png')

end