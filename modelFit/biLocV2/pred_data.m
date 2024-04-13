function [mean_conf, std_mean_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = pred_data(best_para, model_slc, num_rep, sub_slc, ses_slc)

% assign param
aA                    = best_para(1);
bA                    = best_para(2);
sigV1                 = best_para(3);
sigA                  = best_para(4);
sigV2                 = best_para(5);
sigVs                 = [sigV1, sigV2];
sigP                  = best_para(6);
pCommon               = best_para(7);
sigM                  = best_para(8);
cA                    = best_para(9);
cV                    = best_para(10);

% expt info
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
% raw_diff                 = zeros(length(sA), length(sV));
% for ii                 = 1:length(sA)
%     for j                 = 1:length(sV)
%         raw_diff(ii, j)           = sA(ii) - sV(j);
%     end
% end
% abs_diff             = unique(abs(raw_diff))';

lapse                 = 0.02;
muP                   = 0;

[pred_loc, pred_conf] = deal(NaN(num_s, num_cue, num_rel, num_rep));

for j                 = 1:num_s
    for v                 = 1:numel(sigVs)
        [pred_loc(j,:,v,:), pred_conf(j,:,v,:)] = sim_loc_conf(num_rep, sAV(1,j), sAV(2,j),...
            aA, bA, sigA, sigVs(v), muP, sigP, pCommon, sigM, cA, cV, lapse, model_slc);
    end
end

org_pred_loc               = reshape(pred_loc, [numel(sA), numel(sV), num_cue, num_rel, num_rep]);
org_pred_conf              = reshape(pred_conf, [numel(sA), numel(sV), num_cue, num_rel, num_rep]);

% checkFakeData=1;
% 
% if checkFakeData
% 
%     uni_loc = zeros(size(org_pred_loc));
% 
%     loc_a = repmat(sA',[1,numel(sV)]);
%     loc_v = repmat(sV,[numel(sA),1]);
% 
%     uni_loc(:,:,1,1,:) = repmat(loc_a, [1, 1, 1, 1, num_rep]);
%     uni_loc(:,:,2,1,:) = repmat(loc_v, [1, 1, 1, 1, num_rep]);
% 
%     uni_loc(:,:,1,2,:) = uni_loc(:,:,1,1,:);
%     uni_loc(:,:,2,2,:) = uni_loc(:,:,2,1,:);
% 
%     % loc at uni minus loc at bi
%     ve =  mean(uni_loc, 5) - mean(org_pred_loc,5);
% 
%     % diff x cue x reliability
%     [ve_by_raw_diff, all_raw_diffs] = org_by_raw_diffs_4D(ve, sA);
% 
%     % assume participants localized perfectly in the unisensory
%     % condition
%     figure; hold on
%     t = tiledlayout(2, 1);
%     title(t,sprintf('M%s, rep: %i', model_slc, num_rep))
%     xlabel(t, 'Audiovisual discrepancy (deg)');
%     ylabel(t, 'Shift of localization');
%     t.TileSpacing = 'compact';
%     t.Padding = 'compact';
% 
%     for cue = 1:num_cue
%         nexttile
%         title(cue_label{cue})
%         axis equal
%         hold on
% 
%         for rel = 1: numel(sigVs)
% 
%             i_ve = squeeze(ve_by_raw_diff(:, cue, rel));
%             plot(all_raw_diffs, i_ve)
% 
%         end
%         xticks(all_raw_diffs)
%     end
% 
% end

%% load real data

[~, ~, ~, ExpInfo, ~, ScreenInfo] = org_data(sub_slc,ses_slc,'biLoc');
[uni_resp, uni_conf] = org_data(sub_slc,[],'uniLoc');

deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
remapped_sV = ExpInfo.targetPixel .* deg_per_px;

%% analyze predicted bimodal confidence

% percentage of reporting confidence in unimodal task
uni_pconf = mean(mean(uni_conf,3),2); % condition (A,V1,V2) x loc (4)

% organize bi-confidence by discrepancy
[conf_by_diff, abs_diff] = org_by_diffs(org_pred_conf, sA); % {diff} cue x reliability x rep

[mean_conf, std_mean_conf] = deal(nan(numel(abs_diff), 2, 2));
% calculate mean and sd
for cue = 1:2
    for rel = 1:2
        for diff = 1:numel(abs_diff)
            est = squeeze(conf_by_diff{diff}(cue, rel, :))';
            p = sum(est)/numel(est);
            mean_conf(diff, cue, rel) = p;
            std_mean_conf(diff, cue, rel) = sqrt((p*(1-p))/numel(est));
        end
    end
end

%% analyze predicted ventriloquist effect

mean_uni_resp = mean(uni_resp, 3); % condition (A,V1,V2) x loc (4) x rep
loc_a = repmat(mean_uni_resp(1,:)', [1, numel(remapped_sV)]);
remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);

% reshape into dimensions of bi
uni_loc(:,:,1,1:2) = repmat(loc_a, [1,1,1,2]);
uni_loc(:,:,2,1:2) = repmat(remap_loc_v, [1,1,1,2]);

% loc at uni minus loc at bi
ve = uni_loc - mean(org_pred_loc, 5);
std_ve = std(org_pred_loc,[], 5);

% diff x cue x reliability
[mean_ve, raw_diff] = org_by_raw_diffs_4D(ve, sA);
[std_ve, ~] = org_by_raw_diffs_4D(std_ve, sA);

end