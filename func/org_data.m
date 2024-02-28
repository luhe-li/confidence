function [org_resp, org_conf, org_err, ExpInfo] = org_data(sub_slc,ses_slc,exp)
% exp can be 'biLoc', 'uniLoc', or 'pointTask'

[org_resp, org_conf, org_err] = deal([]);

for i = 1:numel(ses_slc)
    % manage path
    cur_dir      = pwd;
    [project_dir, ~]= fileparts(fileparts(cur_dir));
    data_dir     = fullfile(project_dir, 'data',exp);

    % organize data
    switch exp
        case 'biLoc'
            flnm        = sprintf('biLoc_sub%i_ses%i.mat', sub_slc, ses_slc(i));
        case 'uniLoc'
            flnm        = sprintf(['uniLoc_sub%i_ses-' ses_slc], sub_slc);
        case 'pointTask'
            flnm        = sprintf('uniLoc_sub%i_ses-Pointing', sub_slc);
    end
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
   
    % conditions
    seq         = ExpInfo.randAVIdx;
    audIdx      = ExpInfo.audIdx;
    visIdx      = ExpInfo.visIdx;
    cueIdx      = ExpInfo.cueIdx;
    visReliIdx  = ExpInfo.visReliIdx;
    num_rep     = ExpInfo.nRep;
    
    % data
    target      = [Resp.target_cm] * ScreenInfo.numPixels_perCM; % convert target to pixel, center as 0
    resp        = [Resp.response_pixel] - ScreenInfo.xmid; % rescale response with center as 0
    err         = abs(resp - target);
    conf        = [Resp.conf];
    
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