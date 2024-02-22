function [org_resp, org_conf, org_err,ExpInfo] = org_data(sub_slc,ses_slc,total_num_rep)
[org_resp, org_conf, org_err] = deal([]);
for i = 1:numel(ses_slc)
    % manage path
    cur_dir      = pwd;
    [project_dir, ~]= fileparts(fileparts(cur_dir));
    % out_dir      = fullfile(cur_dir, 's1Fig');
    data_dir     = fullfile(project_dir, 'confidence','data','biLoc');
    % if ~exist(out_dir,'dir') mkdir(out_dir); end

    % organize data
    flnm        = sprintf('biLoc_sub%i_ses%i', sub_slc, ses_slc(i));
    load(fullfile(data_dir, flnm))

    % experiment info
    aud_locs    = ExpInfo.speakerLocPixel(ExpInfo.audIdx);
    vis_locs    = round(ExpInfo.targetPixel);
    diffs       = zeros(length(aud_locs), length(vis_locs));
    for i = 1:length(aud_locs)
        for j = 1:length(vis_locs)
            diffs(i, j) = aud_locs(i) - aud_locs(j);
        end
    end
    disc_locs   = unique(abs(diffs));

    % conditions
    seq         = ExpInfo.randAVIdx;
    audIdx      = ExpInfo.audIdx;
    visIdx      = ExpInfo.visIdx;
    cueIdx      = ExpInfo.cueIdx;
    visReliIdx  = ExpInfo.visReliIdx;
    num_rep     = ExpInfo.nRep;
    disc        = abs(seq(1,:) - seq(2,:));
    discIdx     = unique(disc);
    cue_label   = {'Post-cue: A','Post-cue: V'};
    rel_label   = {'High visual reliability','Low visual reliability'};

    % data
    target      = [Resp.target_cm] * ScreenInfo.numPixels_perCM; % convert target to pixel, center as 0
    resp        = [Resp.response_pixel] - ScreenInfo.xmid; % rescale response with center as 0
    err         = abs(resp - target);
    conf        = [Resp.conf_radius_cm];

    % organize response
    [org_resp_temp, org_conf_temp, org_err_temp] = deal(NaN(numel(audIdx), numel(visIdx), numel(cueIdx), numel(visReliIdx), num_rep));
    for trial = 1:length(resp)
        % Find indices for each condition
        aIdx = find(audIdx == seq(1, trial));
        vIdx = find(visIdx == seq(2, trial));
        cIdx = find(cueIdx == seq(3, trial));
        rIdx = find(visReliIdx == seq(4, trial));

        % Find the repetition index for the current combination
        repMatrix = squeeze(org_resp_temp(aIdx, vIdx, cIdx, rIdx, :)); % Squeeze to remove singleton dimensions
        repIdx = find(isnan(repMatrix), 1); % Find the first NaN value

        % Insert the response
        org_resp_temp(aIdx, vIdx, cIdx, rIdx, repIdx) = resp(trial);
        org_conf_temp(aIdx, vIdx, cIdx, rIdx, repIdx) = conf(trial);
        org_err_temp(aIdx, vIdx, cIdx, rIdx, repIdx) = err(trial);
    end
    org_resp = cat(5,org_resp,org_resp_temp);
    org_conf = cat(5,org_conf,org_conf_temp);
    org_err = cat(5,org_err,org_err_temp);
end

end