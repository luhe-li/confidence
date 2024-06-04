function [mean_conf, std_mean_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc, pred)

% mean_conf, std_mean_conf: abs_diff x cue x reliability
% uni_pconf: condition (A,V1,V2) x abs_diff
% mean_ve, std_ve: raw_diff x cue x reliability

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
addpath(genpath(fullfile(project_dir,'func')))

% organize data
if ~exist('pred' ,'var')
    [bi_resp, bi_conf, ~, ExpInfo, ~, ScreenInfo] = org_data(sub_slc,ses_slc,'biLoc');
    audIdx  = ExpInfo.audIdx;
    targetPixel = ExpInfo.targetPixel;
else 
    bi_resp = pred.bi_resp;
    bi_conf = pred.bi_conf;
    audIdx  = pred.audIdx;
    targetPixel = pred.targetPixel;
end    
[uni_resp, uni_conf, ~, uniExpInfo] = org_data(sub_slc,[],'uniLoc');

%% reorganize uni and bi confidence

% sV for bimodal is remapped. Convert them from pixel to VA
% deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
deg_per_px = uniExpInfo.LRmostVisualAngle * 2 / 1024;
sA    = uniExpInfo.speakerLocVA(audIdx);
remapped_sV = targetPixel .* deg_per_px;

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
            std_mean_conf(diff, cue, rel) = std(est)/sqrt(numel(est));
        end
    end
end

%% reorganize ventriloquist effect

% reorganize audlevel for 11 and 12 because they tested 6 uni locations
if sub_slc == 11 || sub_slc == 12
    a_loc_slc = 2:5;
else
    a_loc_slc = 1:4;
end

% if tested auditory locations are the same between uni and bi modal,
% subtract uni_response for VE
if  sum(audIdx == uniExpInfo.audLevel(a_loc_slc)) == 4 % this condition does not interpolate auditory responses
    mean_uni_resp = mean(uni_resp(:,a_loc_slc,:), 3); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(mean_uni_resp(1,:)', [1, numel(remapped_sV)]);
    remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);
% if different, they are not comparable, use interpolated unimodal
% locations
else  % this condition interpolates auditory responses
    % load([sprintf('AVbias_sub%i', ExpInfo.subjID) '.mat'])
    % coefsA = squeeze(Transfer.degCoeff(1, :));
    % There is no need to load bias data because we can calculate it here.
    % This allows fake unimodal data to be used.
    mean_uni_resp = mean(uni_resp, 3);
    sArep = repmat(uniExpInfo.speakerLocVA(uniExpInfo.audLevel)',1,uniExpInfo.nRep);
    uni_a_resp = squeeze(uni_resp(1,:,:));
    mdlA = fitlm(sArep(:),uni_a_resp(:));
    coefsA = table2array(mdlA.Coefficients(:,1));
    interpolate_a_resp = uniExpInfo.speakerLocVA(audIdx) .* coefsA(2) + coefsA(1); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(interpolate_a_resp',[1, numel(remapped_sV)]);
    remap_loc_v = repmat(remapped_sV, [numel(sA), 1]);
    % useful line for debug bias:
    % figure
    % plot(uniExpInfo.speakerLocVA(uniExpInfo.audLevel),mean_uni_resp(1,:),'-o'); hold on; plot(uniExpInfo.speakerLocVA(audIdx),interpolate_a_resp,'-o')
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
