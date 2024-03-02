function sigma_m = get_point_sigM(sub_slc)
% Measure Motor Error from Pointing Task

%% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
data_dir     = fullfile(project_dir, 'data', 'pointTask');
flnm        = sprintf('point_sub%i_ses-Pointing', sub_slc);
load(fullfile(data_dir, flnm));
%%
err = [sortedResp.target_pixel] - [sortedResp.response_pixel] + ScreenInfo.xmid;
sigma_m = std(err);
end
