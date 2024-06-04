function sigma_m = get_point_sigM(sub_slc)
% Measure Motor Error from Pointing Task

%% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
data_dir     = fullfile(project_dir, 'data', 'pointTask');
flnm        = sprintf('point_sub%i_ses-Pointing', sub_slc);
load(fullfile(data_dir, flnm));

%% Calculate motor noise
err = [sortedResp.target_pixel] - [sortedResp.response_pixel] + ScreenInfo.xmid;

% First 4 participants ran into bug that the first trial was skipped, hence
% the response of the first trial was omitted
if sub_slc <= 4
    err = err(2:end);
end
deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
err_deg = err .* deg_per_px;
sigma_m = sqrt(sum((err_deg - mean(err_deg)).^2/(numel(err_deg))));

end
