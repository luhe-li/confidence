
clear; close all; rng('shuffle');

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
result_dir            = fullfile(cur_dir, 'modelFit');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(result_dir))

%% plot set up

clt = [30, 120, 180; % blue
    162, 62, 72;  % dark red
    227, 135, 158]./255; % light red

lw = 2;
datalw = 4;
ftsz = 20;
titlesz = 25;

%% free param

% free parameters
aA = 1;
bA = -0.5;
sigV1 = 1;
sigA = 5;
sigV2 = 4;
sigP = 20;
pCommon = 0.75;
sigC = 2;

%% model info

sA = -10:1:10;
sV = sA;

muP = 0;
lapse = 0.02;
sigVs = [sigV1, sigV2];
num_rep = 1000;
fixP.sigMotor = 1.36;
fixP.bi_nrep = num_rep;
models = {'Heuristic','Suboptimal','Optimal'};
num_models = numel(models);
cue_label = {'Post-cue: A','Post-cue: V'};
num_cue = numel(cue_label);
rel_label = {'High visual reliability','Low visual reliability'};

for d = 1:3

    fixP.model_ind = d;

    % set criterion based on other free parameters
    M_fixP.sA = sA;
    M_fixP.sV = sV;
    M_fixP.model_ind = d;
    M_fixP.num_rep = num_rep;

    [lb, ub] = findConfRange(aA, bA, sigA, sigV1, sigV2, sigP, pCommon, M_fixP);
    intervals = linspace(lb, ub, 5);
    c1 = intervals(2);
    c2 = intervals(3);
    c3 = intervals(4);

    for a = 1:numel(sA)

        for v = 1:numel(sV)

            for r = 1:numel(sigVs)

                fixP.bi_sA = sA(a);
                fixP.bi_sV = sV(v);

                [org_loc{d}(a,v,:,r,:), org_conf{d}(a,v,:,r,:)] = simAllModels(...
                    aA, bA, sigA, sigVs(r), muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP);

            end

        end

    end

end

%% plot

pos = {[0.04 0.08 0.4 0.8],[0.5 0.08 0.4 0.8]};

for d = 1:3

    % organize confidence data: {diff} cue x reliability x rep
    [conf_by_diff, abs_diff] = org_by_diffs(org_conf{d}, sA);

    figure;

    for cue = 1:num_cue

        subplot(1,2,cue)
        set(gca,'Position',pos{cue}, 'LineWidth', lw, 'FontSize', ftsz,'TickDir', 'out')
        set(gcf, 'Position',[0 0 800 350])
        hold on
        title(cue_label{cue},'FontSize',titlesz)
        xticks(round(abs_diff))
        ylim([1 4])
        yticks(1:4)
        xticks(0:5:20)
        xlim([min(abs_diff), max(abs_diff)])

        for rel = 1:numel(sigVs)
            [p_conf, se_conf] = deal(NaN(1, numel(abs_diff)));
            for diff = 1:numel(abs_diff)
                i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
                p = sum(i_conf)/(numel(i_conf));
                p_conf(diff) = p;
                se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
            end
            plot(abs_diff, p_conf, 'Color',clt(rel+1,:),'LineWidth',datalw);
        end

        if cue == 1
            ylabel('Confidence rating','FontSize',titlesz);
            xlabel('Absolute audiovisual discrepancy (deg)','FontSize',titlesz);
        end


    end

        saveas(gca, fullfile(out_dir, sprintf('sim_M%i', d)), 'pdf')

end