function [org_resp, org_conf, org_err, ExpInfo, org_sigVs, ScreenInfo] = org_data(sub_slc,ses_slc,exp)

[org_resp, org_conf, org_err, org_sigVs] = deal([]);

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
            deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
            % data
            target      = [Resp.target_deg]; % convert target to pixel, center as 0
            resp        = [Resp.response_deg]; % rescale response with center as 0
            err         = resp - target; % positive = right; negative = left;
            conf        = [Resp.conf];
            visDotsCoords = NaN(length(Resp),2,10);
            for i = 1:length(Resp)
                visDotsCoords(i,:,:) = Resp(i).vStimDotsCoor;
            end
            sigVs = std(squeeze(visDotsCoords(:,1,:)),[],2) .* deg_per_px;
            % organize response
            [org_resp_temp, org_conf_temp, org_err_temp, org_sigVs_temp] = deal(NaN(numel(audIdx), numel(visIdx), numel(cueIdx), numel(visReliIdx), num_rep));
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
                org_sigVs_temp(aIdx, vIdx, cIdx, rIdx, repIdx) = sigVs(trial);
            end
            org_resp = cat(5,org_resp,org_resp_temp);
            org_conf = cat(5,org_conf,org_conf_temp);
            org_err = cat(5,org_err,org_err_temp);
            org_sigVs = cat(5,org_sigVs,org_sigVs_temp);
            
        end

    case 'uniLoc'

        files = dir(data_dir);
        selectedFiles = files(contains({files.name}, sprintf('uniLoc_sub%i_', sub_slc)));

        % A
        load(fullfile(data_dir, selectedFiles(1).name));
        deg_per_CM  = ExpInfo.LRmostVisualAngle / ExpInfo.LRmostSpeakers2center;
        a_org_target = reshape([sortedResp.target_cm] .* deg_per_CM,[ExpInfo.nRep, ExpInfo.nLevel])'; % convert target to pixel, center as 0
        temp = reshape([sortedResp.response_deg],[ExpInfo.nRep, ExpInfo.nLevel])';
        org_resp(1,:,:) = temp;
        org_err(1,:,:) = temp - a_org_target;
        org_conf(1,:,:) = reshape([sortedResp.conf],[ExpInfo.nRep, ExpInfo.nLevel])';

        % V
        load(fullfile(data_dir, selectedFiles(2).name));
        v1_org_target = reshape([sortedReli1Resp.target_deg],[ExpInfo.nRep, ExpInfo.nLevel])';
        temp2 = reshape([sortedReli1Resp.response_deg],[ExpInfo.nRep, ExpInfo.nLevel])' ;
        org_resp(2,:,:) = temp2;
        org_err(2,:,:) = temp2 - v1_org_target;
        org_conf(2,:,:) = reshape([sortedReli1Resp.conf],[ExpInfo.nRep, ExpInfo.nLevel])';

        v2_org_target = reshape([sortedReli2Resp.target_deg],[ExpInfo.nRep, ExpInfo.nLevel])';
        temp3 = reshape([sortedReli2Resp.response_deg],[ExpInfo.nRep, ExpInfo.nLevel])';
        org_resp(3,:,:) = temp3;
        org_err(3,:,:) = temp3 - v2_org_target;
        org_conf(3,:,:) = reshape([sortedReli2Resp.conf],[ExpInfo.nRep, ExpInfo.nLevel])';

    case 'pointTask'
        flnm        = sprintf('point_sub%i_ses-Pointing', sub_slc);
        load(fullfile(data_dir, flnm));
end


end