
clear; clc; close all; rng('Shuffle');

sub = 4;
ses = 1:2;
models = {'Heuristic','Suboptimal','Optimal'};
m = 3;

flnm = sprintf('FitResults_%s_sub%i_ses%i-%i', models{m}, sub, min(ses), max(ses));

%% manage path
out_dir               = fullfile(pwd, mfilename);

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;
%     251, 154, 153]./255; % light red

%% load best fitting parameters

% look for all the files, load the most updated one that match the selected name
files = dir(fullfile(pwd, [flnm '*.mat'])); % Lists all .mat files in the current folder
load(files(end).name);

% use the best-fitting parameter, simulate responses
p = model.bestP;
num_para = numel(p);

% assign param
aA                    = p(1);
bA                    = p(2);
sigV1                 = p(3);
sigA                  = p(4);
sigV2                 = p(5);
sigVs                 = [sigV1, sigV2];
sigP                  = p(6);
pCommon               = p(7);
sigM                  = p(8);
cA                    = p(9);
cV                    = p(10);

%% simulate prediction of confidence and localization

% expt info
num_rep = model.num_rep;
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};
num_rel               = numel(rel_label);

aud_level             = [6 8 9 11];
aud_VA                = -30:4:30;
sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA;
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of conditions
raw_diff                 = zeros(length(sA), length(sV));
for ii                 = 1:length(sA)
    for j                 = 1:length(sV)
        raw_diff(ii, j)           = sA(ii) - sV(j);
    end
end
abs_diff             = unique(abs(raw_diff))';


lapse                 = 0.02;
muP                   = 0;

[pred_loc, pred_conf] = deal(NaN(num_s, num_cue, num_rel, num_rep));
% num_s: stimulus location combination
% 3 confidence decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep

for j                 = 1:num_s
    for v                 = 1:numel(sigVs)
        [pred_loc(j,:,v,:), pred_conf(j,:,v,:)] = sim_loc_conf(num_rep, sAV(1,j), sAV(2,j),...
            aA, bA, sigA, sigVs(v), muP, sigP, pCommon, sigM, cA, cV, lapse, m);
    end
end

%% organize prediction and data

org_pred_loc               = reshape(pred_loc, [numel(sA), numel(sV), num_cue, num_rel, num_rep]);
org_pred_conf              = reshape(pred_conf, [numel(sA), numel(sV), num_cue, num_rel, num_rep]);

org_loc = data.org_resp;
org_conf = data.org_conf;

%% check localization pred

% load uniloc data
uni_loc = zeros(size(org_pred_loc));

loc_a = repmat(sA',[1,numel(sV)]);
loc_v = repmat(sV,[numel(sA),1]);

uni_loc(:,:,1,1,:) = repmat(loc_a, [1, 1, 1, 1, num_rep]);
uni_loc(:,:,2,1,:) = repmat(loc_v, [1, 1, 1, 1, num_rep]);

uni_loc(:,:,1,2,:) = uni_loc(:,:,1,1,:);
uni_loc(:,:,2,2,:) = uni_loc(:,:,2,1,:);

% loc at uni minus loc at bi
ve =  mean(uni_loc, 5) - mean(org_pred_loc,5);

% diff x cue x reliability
[pred_ve_by_raw_diff, all_raw_diffs] = org_by_raw_diffs_4D(ve, sA);

% assume participants localized perfectly in the unisensory
% condition
figure; hold on
t = tiledlayout(2, 1);
title(t,sprintf('%s, rep: %i', models{m}, num_rep))
xlabel(t, 'Audiovisual discrepancy (deg)');
ylabel(t, 'Shift of localization');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:num_cue
    nexttile
    title(cue_label{cue})
    axis equal
    hold on

    for rel = 1:num_rel

        i_ve = squeeze(pred_ve_by_raw_diff(:, cue, rel));
        plot(all_raw_diffs, i_ve, 'Color',clt(rel+1,:))

    end
    xticks(all_raw_diffs)
end

%% check confidence pred

% organize localization data: {diff} cue x reliability x rep
[pred_conf_by_diff, all_diffs] = org_by_diffs(org_pred_conf, sA);

figure; hold on
t = tiledlayout(2, 1);
title(t,sprintf('%s, rep: %i', models{m}, num_rep))
xlabel(t, 'Absolute audiovisual discrepancy (deg)');
ylabel(t, 'Proportion of confidence report');
t.TileSpacing = 'compact';
t.Padding = 'compact';

for cue = 1:num_cue
    nexttile
    title(cue_label{cue})
    hold on
    for rel = 1:num_rel
        [p_conf, se_conf] = deal(NaN(1, numel(abs_diff)));
        for diff = 1:numel(abs_diff)
            i_conf = squeeze(pred_conf_by_diff{diff}(cue, rel, :))';
            p = sum(i_conf)/numel(i_conf);
            p_conf(diff) = p;
            se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
        end
        plot(abs_diff, p_conf, 'Color',clt(rel+1,:));
        ylim([-0.1, 1.1])
    end
    xticks(abs_diff)
end