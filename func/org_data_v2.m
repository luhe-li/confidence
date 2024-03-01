function [org_resp, org_conf, org_err, ExpInfo] = org_data(sub_slc,ses_slc,exp)

[org_resp, org_conf, org_err] = deal([]);

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
data_dir     = fullfile(project_dir, 'data', exp);

% organize data
switch exp

    case 'biLoc'

        for i = 1:numel(ses_slc)

            flnm        = sprintf('biLoc_sub%i_ses%i.mat', sub_slc, ses_slc(i));
            load(fullfile(data_dir, flnm));

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
            err         = resp - target; % positive = right; negative = left;
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

    case 'uniLoc'

        files = dir(data_dir);
        selectedFiles = files(contains({files.name}, sprintf('uniLoc_sub%i_', sub_slc)));

        % A
        load(fullfile(data_dir, selectedFiles(1).name));

        a_org_target = reshape([sortedResp.target_cm],[ExpInfo.nRep, ExpInfo.nLevel]) * ScreenInfo.numPixels_perCM; % convert target to pixel, center as 0
        
        temp = reshape([sortedResp.response_pixel],[ExpInfo.nRep, ExpInfo.nLevel]) - ScreenInfo.xmid;
        org_resp(1,:,:) = temp;
        org_err(1,:,:) = temp - a_org_target;
        org_conf(1,:,:) = reshape([sortedResp.conf],[ExpInfo.nRep, ExpInfo.nLevel]);

        % V
        load(fullfile(data_dir, selectedFiles(2).name));
        v1_org_target = reshape([sortedReli1Resp.target_pixel],[ExpInfo.nRep, ExpInfo.nLevel]);
        temp2 = reshape([sortedReli1Resp.response_pixel],[ExpInfo.nRep, ExpInfo.nLevel]) - ScreenInfo.xmid;
        org_resp(2,:,:) = temp2;
        org_err(2,:,:) = temp2 - v1_org_target;
        org_conf(2,:,:) = reshape([sortedReli1Resp.conf],[ExpInfo.nRep, ExpInfo.nLevel]);

        v2_org_target = reshape([sortedReli2Resp.target_pixel],[ExpInfo.nRep, ExpInfo.nLevel]);
        temp3 = reshape([sortedReli2Resp.response_pixel],[ExpInfo.nRep, ExpInfo.nLevel]) - ScreenInfo.xmid;
        org_resp(3,:,:) = temp3;
        org_err(3,:,:) = temp3 - v2_org_target;
        org_conf(3,:,:) = reshape([sortedReli2Resp.conf],[ExpInfo.nRep, ExpInfo.nLevel]);

    case 'pointTask'
        flnm        = sprintf('uniLoc_sub%i_ses-Pointing', sub_slc);
end


end