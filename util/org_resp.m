function [org, exp_info] = org_resp(sub_init, all_ses, exp)

% manage path
[project_dir, ~] = fileparts(pwd);
if ~strcmp(fileparts(project_dir), 'confidence')
    [project_dir, ~] = fileparts(project_dir);
end
data_dir     = fullfile(project_dir, 'data', exp);
org = struct();

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
                org.uni_target(j,:) = temp_loc(:, 1);  % Only take the first repetition because locations are the same for all

                % localization estimates data
                temp_est = reshape([sortedResp(1:end).response_cm],[n_rep, n_level])';
                org.uni_loc(j, :, :) = temp_est;

                % Calculate error
                temp_err = temp_est - temp_loc;
                org.uni_err(j, :, :) = temp_err;
                org.uni_abs_err(j, :, :) = abs(temp_err);

                % Calculate mean and standard deviation of estimates
                org.uni_loc_mu(j, :) = mean(temp_est, 2);
                org.uni_loc_sd(j, :) = std(temp_est,[],2);

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

        if ~exist('all_ses', 'var')
            fprintf('You need to specify bimodal sessions.\n');
            return;
        end
        n_ses = length(all_ses);

        % Load data for each session
        for j = 1:n_ses
            session = all_ses(j);
            file_path = fullfile(data_dir, sprintf('biLoc_sub-%s_ses-%d.mat', sub_init, session));
            if exist(file_path, 'file')
                data = load(file_path);
                exp_info = squeeze(data.ExpInfo);
                sortedResp = squeeze(data.sortedResp);

                n_rep = double(exp_info.nRep);
                n_level = double(exp_info.nLevel);

                % Determine the range of trials to use based on the session index
                if j == 1
                    trials_range = 1:floor(length(sortedResp)/2);
                elseif j == 2
                    trials_range = floor(length(sortedResp)/2) + 1:length(sortedResp);
                else
                    error('Unexpected session index: %d', j);
                end

                % target_cm is the true location of the target
                % temp_loc: location x rep
                temp_loc = reshape([sortedResp(trials_range).target_cm], [n_rep, n_level])';
                org.bi_target(j, :) = temp_loc(:, 1);  % Only take the first repetition because locations are the same for all

                % localization estimates data
                temp_est = reshape([sortedResp(trials_range).response_cm], [n_rep, n_level])';
                org.bi_loc(j, :, :) = temp_est;

                % Calculate error
                temp_err = temp_est - temp_loc;
                org.bi_err(j, :, :) = temp_err;
                org.bi_abs_err(j, :, :) = abs(temp_err);

                % Calculate mean and standard deviation of estimates
                org.bi_loc_mu(j, :) = mean(temp_est, 2);
                org.bi_loc_sd(j, :) = std(temp_est, [], 2);

                % Confidence data
                temp_conf = reshape([sortedResp(1:end).conf_radius_cm], [n_rep, n_level])';
                org.bi_conf(j, :, :) = temp_conf;

                % conversion function for plotting
                org.deg_from_cm = @(cm) rad2deg(atan(cm / exp_info.sittingDistance));

            else
                fprintf('File %s not found.\n', file_path);
            end
        end

    case 'pointing'


end


end