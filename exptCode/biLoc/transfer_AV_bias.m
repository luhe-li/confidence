function [coeff_va_slope, coeff_va_intercept,...
    coeff_px_slope, coeff_px_intercept,...
    coeff_cm_slope, coeff_cm_intercept] = transfer_AV_bias(out_dir, subjIni)
%TRANSFER_AV_BIAS Compute or retrieve AV bias coefficients for a subject
%   [coeff_va_slope, coeff_px_slope, coeff_cm_slope] = TRANSFER_AV_BIAS(out_dir, subjIni)
%   Checks if 'AVbias.csv' exists in 'out_dir', if not, creates one.
%   If 'AVbias.csv' exists, loads it. The file should have columns:
%   'subj_init', 'coeff_va_slope', 'coeff_va_intercept',
%   'coeff_px_slope', 'coeff_px_intercept',
%   'coeff_cm_slope', 'coeff_cm_intercept'.
%   If the file does not exist, or if the subject is not in the file,
%   it computes the coefficients.
%   Saves the computed coefficients to the corresponding columns and saves the csv.
%   Outputs the corresponding slopes 'coeff_va_slope', 'coeff_px_slope', 'coeff_cm_slope' for the subject.

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Check if AVbias.csv exists in out_dir
csv_file = fullfile(out_dir, 'AVbias.csv');

if ~exist(csv_file, 'file')
    % File does not exist, create an empty table
    fprintf('AVbias.csv does not exist in %s. Creating a new file.\n', out_dir);
    AVbias = table('Size', [0 7], 'VariableTypes', {'string', 'double', ...
        'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'subj_init', 'coeff_va_slope', ...
        'coeff_va_intercept', 'coeff_px_slope', 'coeff_px_intercept', 'coeff_cm_slope', 'coeff_cm_intercept'});

else
    % Load the existing file
    fprintf('AVbias.csv found in %s. Loading existing data.\n', out_dir);
    AVbias = readtable(csv_file);
end

% Now, check if AVbias has an entry for subjIni
idx = find(strcmp(AVbias.subj_init, subjIni));

if isempty(idx)
    % Entry for this subject does not exist
    fprintf('Subject %s not found in AVbias.csv. Computing new coefficients.\n', subjIni);
    
    % Load data and compute coefficients
    [Resp_dir, ~] = fileparts(fileparts(out_dir));
    addpath(genpath(fullfile(Resp_dir, 'Resp','uniloc')));
    load(sprintf('uniLoc_sub-%s_ses-%s', subjIni,'A'));
    
    % Compute coefficient for degrees
    x_deg = [Resp.target_deg];
    y_deg = [Resp.response_deg];
    coeffs_deg = polyfit(x_deg, y_deg, 1);
    coeff_va_slope = coeffs_deg(1);
    coeff_va_intercept = coeffs_deg(2);

    % Compute coefficient for pixels
    x_px = [Resp.target_pixel];
    y_px = [Resp.response_pixel];
    coeffs_px = polyfit(x_px, y_px, 1);
    coeff_px_slope = coeffs_px(1);
    coeff_px_intercept = coeffs_px(2);

    % Compute coefficient for centimeters
    x_cm = [Resp.target_cm];
    y_cm = [Resp.response_cm];
    coeffs_cm = polyfit(x_cm, y_cm, 1);
    coeff_cm_slope = coeffs_cm(1);
    coeff_cm_intercept = coeffs_cm(2);

    % Save the new coefficients in AVbias table
    new_row = {subjIni, coeff_va_slope, coeff_va_intercept, ...
        coeff_px_slope, coeff_px_intercept, coeff_cm_slope, coeff_cm_intercept};
    AVbias = [AVbias; new_row];
    
    fprintf('Coefficients for subject %s have been computed and added to AVbias.csv.\n', subjIni);
else
    % Entry for this subject exists
    fprintf('Subject %s found in AVbias.csv. Retrieving existing coefficients.\n', subjIni);
    
    % Get the coefficients for this subject
    coeff_va_slope = AVbias.coeff_va_slope(idx);
    coeff_va_intercept = AVbias.coeff_va_intercept(idx);
    coeff_px_slope = AVbias.coeff_px_slope(idx);
    coeff_px_intercept = AVbias.coeff_px_intercept(idx);
    coeff_cm_slope = AVbias.coeff_cm_slope(idx);
    coeff_cm_intercept = AVbias.coeff_cm_intercept(idx);
end

% Save the updated AVbias table to csv
writetable(AVbias, csv_file);
fprintf('AVbias.csv has been updated.\n');

end