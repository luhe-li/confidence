
% understand the nonlinear trend in confidence rating as audiovisual
% discrepancy increases

clear; 

i_gt = [     1,   0,     1,   4,    5,    8,  0.57,   0.5];
ds_loc                = {'Model averaging'};
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

sA                    = -10:10;
sV                    = -10:10;
num_rep = 1000;

clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

for d = 2

    % assign simulation parameters
    num_para              = length(i_gt);
    aA                    = i_gt(1);
    bA                    = i_gt(2);
    sigV1                 = i_gt(3);
    sigA                  = i_gt(4);
    sigV2                 = i_gt(5);
    sigVs                 = [sigV1, sigV2];
    sigP                  = i_gt(6);
    pCommon               = i_gt(7);
    sigC                  = i_gt(8);
    lapse                 = 0.02;
    muP                   = 0;

    %% set criterion based on other free parameters

    M_fixP.sA = sA;
    M_fixP.sV = sV;
    M_fixP.model_ind = d;
    M_fixP.num_rep = num_rep;

    [lb, ub] = findConfRange(aA, bA, sigA, sigV1, sigV2, sigP, pCommon, M_fixP);
    intervals = linspace(lb, ub, 5);
    c1 = intervals(2);
    c2 = intervals(3);
    c3 = intervals(4);

    %% simulate fake data

    % num_s: stimulus location combination
    % 3 confidence decision strategies
    % 2 modalities(1 = aud, 2 = vis)
    % 2 visual reliabilities
    % num_rep
    [org_loc, org_conf] = deal(NaN(numel(sA), numel(sV), num_cue, numel(sigVs), num_rep));

    for a = 1:numel(sA)

        for v = 1:numel(sV)

            for r = 1:numel(sigVs)

                fixP.sA = sA(a);
                fixP.sV = sV(v);
                fixP.model_ind = d;
                fixP.sigMotor = 1.36; % in deg, emperical motor noise averaged from first four participants
                fixP.num_rep = num_rep;

                [org_loc(a,v,:,r,:), org_conf(a,v,:,r,:), variance(a,v,:,r,:), ...
                    norm_var(a,v,:,r,:), est_var(a,v,:,r,:)] = simAllModels(...
                    aA, bA, sigA, sigVs(r), muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP);

            end

        end

    end

end

%% organize by discrepancy

% organize confidence data: {diff} cue x reliability x rep
[var_by_diff, all_diffs] = org_by_diffs(norm_var, sA);

figure; hold on
t = tiledlayout(2, 1);
title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
xlabel(t, 'Absolute audiovisual discrepancy (deg)');
ylabel(t, 'Confidence rating (1-4)');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:num_cue
    nexttile
    title(cue_label{cue})
    hold on
    for rel = 1: numel(sigVs)
        [p_conf, se_conf] = deal(NaN(1, numel(all_diffs)));
        for diff = 1:numel(all_diffs)
            i_conf = squeeze(var_by_diff{diff}(cue, rel, :))';
            p = sum(i_conf)/(numel(i_conf));
            p_conf(diff) = p;
            se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
        end
        plot(all_diffs, p_conf, 'Color',clt(rel+1,:));
        ylim([1,5])
    end
    xticks(all_diffs)
end
