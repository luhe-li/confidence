function [org, exp_info] = org_resp(sub_init, all_ses, exp)

% manage path
[project_dir, ~]= fileparts(pwd);
data_dir     = fullfile(project_dir, 'data', exp);

switch exp
    case 'uniLoc'

        if ~exist('all_ses','var'); all_ses = {'A','V'}; end
        n_ses = length(all_ses);

        % Load data for each session
        for j = 1:n_ses
            session = all_ses{j};
            file_path = fullfile(data_dir, sprintf('uniLoc_sub-%s_ses-%s.mat', sub_init, session));
            if exist(file_path, 'file')
                data = load(file_path);
                exp_info = squeeze(data.ExpInfo);
                sortedResp = squeeze(data.sortedResp);

                n_rep = double(exp_info.nRep);
                n_level = double(exp_info.nLevel);

                % target_cm is the true location of the target
                % temp_loc: location x rep
                temp_loc = reshape([sortedResp(1:end).target_cm],[n_rep, n_level])'; 
                org.uni_target(:) = temp_loc(:, 1);  % Only take the first repetition because locations are the same for all

                % localization estimates data
                temp_est = reshape([sortedResp(1:end).response_cm],[n_rep, n_level])'; 
                org.uni_loc(j, :, :) = temp_est;

                % Calculate error
                temp_err = temp_est - temp_loc;
                org.uni_err(j, :, :) = temp_err;

                % Calculate mean and standard deviation of estimates
                org.uni_loc_mu(j, :) = mean(temp_est, 2);
                org.uni_loc_sd(j, :) = sqrt(sum((temp_est - mean(temp_est,2)).^2)./numel(temp_est));

                % Confidence data
                temp_conf = reshape([sortedResp(1:end).conf_radius_cm],[n_rep, n_level])'; 
                org.uni_conf(j, :, :) = temp_conf;

                % conversion function for plotting
                org.deg_from_cm = @(cm) rad2deg(atan(cm / exp_info.sittingDistance));

            else
                fprintf('File %s not found.\n', file_path);
            end
        end

    case 'biLoc'

    case 'pointing'


end


end