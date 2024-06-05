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

% load unimodal prediction if it exists
if exist('pred.uni_loc','var')
    uni_loc = pred.uni_loc;
    uni_conf = pred.uni_conf;
else % load unimodal data if there's no prediction
    [uni_loc, uni_conf, ~, ~] = org_data(sub_slc,[],'uniLoc');
    if exist('pred','var') % load bimodal prediction if it exists
        bi_loc = pred.bi_loc;
        bi_conf = pred.bi_conf;
        biExpInfo = pred.biExpInfo;
    else % load bimodal data if there's no prediction
        [bi_loc, bi_conf, ~, biExpInfo] = org_data(sub_slc,ses_slc,'biLoc');
        [uni_loc, uni_conf, ~, ~] = org_data(sub_slc,[],'uniLoc');
    end
end

% load uni data anyway

% get stimulus loactions
% get stimulus loactions;
bi_sA = unique(biExpInfo.randAudVA);
bi_sV = unique(biExpInfo.randVisVA);
uni_sA = bi_sA;

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

% interpolates auditory responses as a straight line
if interpolateUni || ~isequal(uni_sA, bi_sA)

    % load([sprintf('AVbias_sub%i', ExpInfo.subjID) '.mat'])
    % coefsA = squeeze(Transfer.degCoeff(1, :));
    % There is no need to load bias data because we can calculate it here.
    % This allows fake unimodal data to be used.
    mean_uni_resp = mean(uni_loc, 3);
%     sArep = repmat(bi_sA',1,uniExpInfo.nRep);
    sArep = repmat(bi_sA',1,20); 
    uni_a_resp = squeeze(uni_loc(1,:,:));
    mdlA = fitlm(sArep(:),uni_a_resp(:));
    coefsA = table2array(mdlA.Coefficients(:,1));
    interpolate_a_resp = uni_sA .* coefsA(2) + coefsA(1); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(interpolate_a_resp',[1, numel(bi_sV)]);
    loc_v = repmat(bi_sV, [numel(bi_sA), 1]);
    % useful line for debug bias:
    % figure
    % plot(uniExpInfo.speakerLocVA(uniExpInfo.audLevel),mean_uni_resp(1,:),'-o'); hold on; plot(uniExpInfo.speakerLocVA(audIdx),interpolate_a_resp,'-o')

% use mean uni loc responses 
else

    mean_uni_resp = mean(uni_loc, 3); % condition (A,V1,V2) x loc (4)
    loc_a = repmat(mean_uni_resp(1,:)', [1, numel(bi_sV)]);
    loc_v = repmat(bi_sV, [numel(bi_sA), 1]);

end

% reshape into dimensions of bi
org_uni_loc(:,:,1,1:2) = repmat(loc_a, [1,1,1,2]);
org_uni_loc(:,:,2,1:2) = repmat(loc_v, [1,1,1,2]);

% loc at uni minus loc at bi
ve = mean(bi_loc, 5) - org_uni_loc;
std_ve = std(bi_loc,[], 5) ./ sqrt(size(bi_loc,5));

% diff x cue x reliability
[mean_ve, raw_diff] = org_by_raw_diffs_4D(ve, bi_sA);
[std_ve, ~] = org_by_raw_diffs_4D(std_ve, bi_sA);

end
