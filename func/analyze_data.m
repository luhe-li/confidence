function [mean_conf, std_mean_conf, uni_pconf, abs_diff,...
    mean_ve, std_ve, raw_diff] = analyze_data(sub_slc, ses_slc, pred, interpolateUni)

% mean_conf, std_mean_conf: abs_diff x cue x reliability
% uni_pconf: condition (A,V1,V2) x abs_diff
% mean_ve, std_ve: raw_diff x cue x reliability

if ~exist('interpolateUni' ,'var'); interpolateUni = true; end

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
addpath(genpath(fullfile(project_dir,'func')))

% load bimodal data if it doesn't exist
if ~exist('pred' ,'var')
    [bi_loc, bi_conf, ~, biExpInfo] = org_data(sub_slc,ses_slc,'biLoc');

% or use predicted data
else 
    bi_loc = pred.bi_loc;
    bi_conf = pred.bi_conf;
    biExpInfo = pred.biExpInfo;
end    
% load uni data anyway
[uni_resp, uni_conf, ~, uniExpInfo] = org_data(sub_slc,[],'uniLoc');

% get stimulus loactions
uni_sA = unique(uniExpInfo.randAudVA);
uni_sV = unique(uniExpInfo.randVisVA);
bi_sA = unique(biExpInfo.randAudVA);
bi_sV = unique(biExpInfo.randVisVA);

%% reorganize uni and bi confidence

% percentage of reporting confidence in unimodal task
uni_pconf = mean(mean(uni_conf,3),2); % condition (A,V1,V2) x loc (4)

% organize bi-confidence by discrepancy
[conf_by_diff, abs_diff] = org_by_diffs(bi_conf, bi_sA); % {diff} cue x reliability x rep

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
    bi_idx = 2:5;
else
    bi_idx = 1:4;
end

% interpolates auditory responses as a straight line
if interpolateUni || ~isequal(uni_sA, bi_sA)

    % load([sprintf('AVbias_sub%i', ExpInfo.subjID) '.mat'])
    % coefsA = squeeze(Transfer.degCoeff(1, :));
    % There is no need to load bias data because we can calculate it here.
    % This allows fake unimodal data to be used.
    mean_uni_resp = mean(uni_resp, 3);
    sArep = repmat(uniExpInfo.speakerLocVA(uniExpInfo.audLevel)',1,uniExpInfo.nRep);
    uni_a_resp = squeeze(uni_resp(1,:,:));
    mdlA = fitlm(sArep(:),uni_a_resp(:));
    coefsA = table2array(mdlA.Coefficients(:,1));
    interpolate_a_resp = uni_sA(bi_idx) .* coefsA(2) + coefsA(1); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(interpolate_a_resp',[1, numel(bi_sV)]);
    loc_v = repmat(bi_sV, [numel(bi_sA), 1]);
    % useful line for debug bias:
    % figure
    % plot(uniExpInfo.speakerLocVA(uniExpInfo.audLevel),mean_uni_resp(1,:),'-o'); hold on; plot(uniExpInfo.speakerLocVA(audIdx),interpolate_a_resp,'-o')

% use mean uni loc responses 
else

    mean_uni_resp = mean(uni_resp(:,bi_idx,:), 3); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(mean_uni_resp(1,:)', [1, numel(bi_sV)]);
    loc_v = repmat(bi_sV, [numel(bi_sA), 1]);

end

% reshape into dimensions of bi
uni_loc(:,:,1,1:2) = repmat(loc_a, [1,1,1,2]);
uni_loc(:,:,2,1:2) = repmat(loc_v, [1,1,1,2]);

% loc at uni minus loc at bi
ve = mean(bi_loc, 5) - uni_loc;
std_ve = std(bi_loc,[], 5);

% diff x cue x reliability
[mean_ve, raw_diff] = org_by_raw_diffs_4D(ve, bi_sA);
[std_ve, ~] = org_by_raw_diffs_4D(std_ve, bi_sA);

end
