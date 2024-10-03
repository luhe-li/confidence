function [opt_est, opt_radius, opt_gain] = eGain_optimized(myPDF, maxScore, minScore, elbow, center_axis)
% Optimized version of the eGain function

% Normalize myPDF along the second dimension
myPDF = myPDF ./ sum(myPDF, 2);
[n_trial, n_est] = size(myPDF);
step = center_axis(2) - center_axis(1);
idx_elbow = round(elbow / step);
penaltyRange = maxScore - minScore;

% Precompute cumulative sums
cumPDF = cumsum(myPDF, 2); % size n_trial x n_est
cumPDF_ext = [zeros(n_trial, 1), cumPDF]; % size n_trial x (n_est + 1)

% Initialize variables
maxGain = NaN(n_trial, n_est);
optRadius = NaN(n_trial, n_est);

% Loop over possible location estimates
for idx_est = 1:n_est
    % Maximum possible confidence radius given the estimate position
    idx_conf = min([idx_est - 1, n_est - idx_est]);
    
    if idx_conf < 1
        maxGain(:, idx_est) = 0;
        optRadius(:, idx_est) = 0;
    else
        % Possible confidence radii
        confRadius = 0:idx_conf;
        lengthRatio = confRadius / idx_elbow;
        costFun = maxScore - lengthRatio * penaltyRange;
        costFun = max(costFun, minScore); % Ensure costFun >= minScore
        
        % Preallocate erCDF and gainFun
        erCDF = zeros(n_trial, length(confRadius));
        
        % Compute erCDF for all trials and radii
        idx_rad_vec = confRadius;
        idx_left = idx_est - idx_rad_vec;
        idx_right = idx_est + idx_rad_vec;
        idx_left = max(1, idx_left);
        idx_right = min(n_est, idx_right);
        idx_left_ext = idx_left;
        idx_right_ext = idx_right + 1; % Adjust for cumPDF_ext indexing
        
        for idx_r = 1:length(idx_rad_vec)
            left_indices = idx_left_ext(idx_r);
            right_indices = idx_right_ext(idx_r);
            if left_indices > 1
                erCDF(:, idx_r) = cumPDF_ext(:, right_indices) - cumPDF_ext(:, left_indices - 1);
            else
                erCDF(:, idx_r) = cumPDF_ext(:, right_indices);
            end
        end
        
        % Compute gainFun for all trials
        gainFun = erCDF .* costFun;
        
        % Find maximum gain and corresponding radius for each trial
        [maxGain(:, idx_est), idx_opt_rad] = max(gainFun, [], 2);
        optRadius(:, idx_est) = confRadius(idx_opt_rad);
    end
end

% Convert optimal radius from axis unit to cm
optRadiusCM = optRadius .* step;

% Find the estimate that gives the maximum gain for each trial
[opt_gain, idx_opt] = max(maxGain, [], 2);
opt_est = center_axis(idx_opt);
opt_radius = arrayfun(@(tt) optRadiusCM(tt, idx_opt(tt)), 1:n_trial).';

end