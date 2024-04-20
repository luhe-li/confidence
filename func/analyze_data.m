function [mean_conf, std_mean_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc)

% mean_conf, std_mean_conf: abs_diff x cue x reliability
% uni_pconf: condition (A,V1,V2) x abs_diff
% mean_ve, std_ve: raw_diff x cue x reliability

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
addpath(genpath(fullfile(project_dir,'func')))

% organize data
[bi_resp, bi_conf, ~, ExpInfo, ~, ScreenInfo] = org_data(sub_slc,ses_slc,'biLoc');
[uni_resp, uni_conf] = org_data(sub_slc,[],'uniLoc');

%% reorganize uni and bi confidence

deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
sA    = ExpInfo.speakerLocVA(ExpInfo.audIdx);
remapped_sV = ExpInfo.targetPixel .* deg_per_px;

% percentage of reporting confidence in unimodal task
uni_pconf = mean(mean(uni_conf,3),2); % condition (A,V1,V2) x loc (4)

% organize bi-confidence by discrepancy
[conf_by_diff, abs_diff] = org_by_diffs(bi_conf, sA); % {diff} cue x reliability x rep

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

%% reorganize ventriloquist effect

mean_uni_resp = mean(uni_resp, 3); % condition (A,V1,V2) x loc (4)
loc_a = repmat(mean_uni_resp(1,:)', [1, numel(remapped_sV)]);
% loc_a = repmat(sA',[1, numel(remapped_sV)]);
remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);

% reshape into dimensions of bi
uni_loc(:,:,1,1:2) = repmat(loc_a, [1,1,1,2]);
uni_loc(:,:,2,1:2) = repmat(remap_loc_v, [1,1,1,2]);

% loc at uni minus loc at bi
ve = mean(bi_resp, 5) - uni_loc;
std_ve = std(bi_resp,[], 5);

% diff x cue x reliability
[mean_ve, raw_diff] = org_by_raw_diffs_4D(ve, sA);
[std_ve, ~] = org_by_raw_diffs_4D(std_ve, sA);

end
