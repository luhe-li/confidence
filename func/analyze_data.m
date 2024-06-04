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
[uni_resp, uni_conf, ~, uniExpInfo] = org_data(sub_slc,[],'uniLoc');

%% reorganize uni and bi confidence

% sV for bimodal is remapped. Convert them from pixel to VA
% deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
deg_per_px = ExpInfo.LRmostVisualAngle * 2 / ScreenInfo.xaxis;
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

% if tested auditory locations are the same between uni and bi modal,
% subtract uni_response for VE
if  sum(ExpInfo.audIdx == uniExpInfo.audLevel) == 4
    mean_uni_resp = mean(uni_resp, 3); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(mean_uni_resp(1,:)', [1, numel(remapped_sV)]);
    remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);
% if different, they are not comparable, use interpolated unimodal
% locations
else
    % load([sprintf('AVbias_sub%i', ExpInfo.subjID) '.mat'])
    % coefsA = squeeze(Transfer.degCoeff(1, :));
    mean_uni_resp = mean(uni_resp, 3);
    sArep = repmat(ExpInfo.speakerLocVA(uniExpInfo.audLevel)',1,uniExpInfo.nRep);
    uni_a_resp = squeeze(uni_resp(1,:,:));
    mdlA = fitlm(sArep(:),uni_a_resp(:));
    coefsA = table2array(mdlA.Coefficients(:,1));
    interpolate_a_resp = ExpInfo.speakerLocVA(ExpInfo.audIdx) .* coefsA(2) + coefsA(1); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(interpolate_a_resp',[1, numel(remapped_sV)]);
    remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);
    % useful line for debug bias:
    % figure
    % plot(ExpInfo.speakerLocVA(uniExpInfo.audLevel),mean_uni_resp(1,:),'-o'); hold on; plot(ExpInfo.speakerLocVA(ExpInfo.audIdx),interpolate_a_resp,'-o')
end

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
