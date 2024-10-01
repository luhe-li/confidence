% plot discrimination data from staircase for each standard location

clear; clc; %close all;

%% set up

sub = 'LL';
ses = 'V';
save_fig = 1;
str = '';

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, mfilename);
addpath(genpath(fullfile(project_dir, 'data',['uniDiscrimination', str])));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load data
load(sprintf('uniDis%s_sub-%s_ses-%s',str, sub, ses));
n_staircase = ExpInfo.n_staircase;

%% Plot the interleaved staircases

% Prepare figure
figure;
hold on;

% Get the number of standard locations
n_standard_locations = size(Resp.comparison_loc, 1);

% Define the number of rows and columns for subplots
n_rows = ceil(sqrt(n_standard_locations));
n_cols = ceil(n_standard_locations / n_rows);

% Loop over each standard location
for k = 1:n_standard_locations
    subplot(n_rows, n_cols, k);
    hold on;
    
    trialNum = 1:ExpInfo.n_trial;
    n_easy_trials = ExpInfo.n_easy_trial_per_s;
    speaker_cm = ExpInfo.speaker_level_cm;
    
    % For staircase == 1
    condition1 = 1; % Assuming condition 1 is in row 1
    comparison_loc1 = Resp.comparison_loc(k, condition1, 1:(end - n_easy_trials));
    correct1 = Resp.correct(k, condition1, 1:(end - n_easy_trials));
    
    % convert comparirson_loc to cm if ses == 'A'
    if strcmp(ses, 'A')
        comparison_values1 = speaker_cm(comparison_loc1);
    else
        comparison_values1 = comparison_loc1;
    end
    
    % Plot comparison locations
    h_line1 = plot(trialNum, comparison_values1, 'b-', 'LineWidth', 1.5);
    
    % Mark correct and incorrect responses
    correctTrials1 = correct1 == 1;
    incorrectTrials1 = correct1 == -1;
    
    h_correct1 = plot(trialNum(correctTrials1), comparison_values1(correctTrials1), 'bo', 'MarkerFaceColor', 'b');
    h_incorrect1 = plot(trialNum(incorrectTrials1), comparison_values1(incorrectTrials1), 'bx', 'LineWidth', 2, 'MarkerSize', 10);
    
    % For staircase == 2
    condition2 = 2; % Assuming condition 2 is in row 2
    comparison_loc2 = Resp.comparison_loc(k, condition2, 1:(end - n_easy_trials));
    correct2 = Resp.correct(k, condition2, 1:(end - n_easy_trials));
    
    % convert comparirson_loc to cm if ses == 'A'
    if strcmp(ses, 'A')
        comparison_values2 = speaker_cm(comparison_loc2);
    else
        comparison_values2 = comparison_loc2;
    end
    
    % Plot comparison locations
    h_line2 = plot(trialNum, comparison_values2, 'r-', 'LineWidth', 1.5);
    
    % Mark correct and incorrect responses
    correctTrials2 = correct2 == 1;
    incorrectTrials2 = correct2 == -1;
    
    h_correct2 = plot(trialNum(correctTrials2), comparison_values2(correctTrials2), 'ro', 'MarkerFaceColor', 'r');
    h_incorrect2 = plot(trialNum(incorrectTrials2), comparison_values2(incorrectTrials2), 'rx', 'LineWidth', 2, 'MarkerSize', 10);
    
    % Add labels and legend
    xlabel('Trial number');
    ylabel('Comparison stimulus location (cm)');
    ylim([-80 80]);
    xlim([0 max(trialNum)]);
    title(sprintf('Standard Location %d', k));
    
    if k == 1
        legend([h_line1, h_line2], '1-up-2-down', '2-up-1-down');
    end
    
    grid on;
    hold off;
end

% Add a main title for the entire figure
sgtitle([ses '-2IFC Discrimination']);
saveas(gcf, fullfile(out_dir, sprintf('stim_check_%s_%s', str, ses)), 'png');

%% plot raw data for PMF
% Determine bin centers and edges
if strcmp(ses, 'A')
    % Convert speaker index to cm
    standard_loc = mean(speaker_cm(ExpInfo.standard_loc));
    comparison_loc = speaker_cm(ExpInfo.comparison_loc);
    Discrepancies = Discrepancies * (speaker_cm(2) - speaker_cm(1));
    bin_width = 5;
else
    standard_loc = ExpInfo.standard_loc;
    comparison_loc = ExpInfo.comparison_loc;
    bin_width = 5;
end

% Determine the maximum absolute discrepancy, rounded up to the nearest multiple of bin_width/2
maxAbsDiscrepancy = ceil(max(abs(comparison_loc)) / (bin_width / 2)) * (bin_width / 2);

% Create bin centers from -maxAbsDiscrepancy to +maxAbsDiscrepancy in steps of bin_width
binCenters = -maxAbsDiscrepancy:bin_width:maxAbsDiscrepancy;

% Compute bin edges by shifting bin centers by half the bin width
binEdges = binCenters - bin_width / 2;
% Ensure that the last edge covers the maximum discrepancy
binEdges = [binEdges, binCenters(end) + bin_width / 2];

% Initialize figure
f1 = figure;
hold on;

% Loop over each standard location
for k = 1:n_standard_locations
    % Bin the discrepancies
    [~, ~, binIdx] = histcounts(Discrepancies, binEdges);

    % Initialize arrays for counts
    numTrialsPerBin = zeros(size(binCenters));
    numRightResponsesPerBin = zeros(size(binCenters));

    % Loop over bins to count responses
    for b = 1:length(binCenters)
        idx = binIdx == b;
        numTrialsPerBin(b) = sum(idx);
        numRightResponsesPerBin(b) = sum(Resp.resp(k, :, idx) == 1); % 1 indicates 'Comparison stimulus is more right'
    end

    % Calculate percentage
    PercentageRight = zeros(size(binCenters));
    PercentageRight(numTrialsPerBin > 0) = (numRightResponsesPerBin(numTrialsPerBin > 0) ./ numTrialsPerBin(numTrialsPerBin > 0)) * 100;

    % Identify bins with at least one trial
    validBins = numTrialsPerBin > 0;

    % Extract valid data
    validBinCenters = binCenters(validBins);
    validPercentageRight = PercentageRight(validBins);
    validNumTrials = numTrialsPerBin(validBins);

    % Scale marker sizes (e.g., map numTrialsPerBin to a reasonable marker size range)
    % Define minimum and maximum marker sizes
    minMarkerSize = 10; % Minimum marker size
    maxMarkerSize = 100; % Maximum marker size

    % Scale numTrialsPerBin to marker sizes between minMarkerSize and maxMarkerSize
    if max(validNumTrials) == min(validNumTrials)
        % All bins have the same number of trials
        markerSizes = ones(size(validNumTrials)) * ((minMarkerSize + maxMarkerSize) / 2);
    else
        % Normalize numTrialsPerBin to range [0,1]
        markerSizes = (validNumTrials - min(validNumTrials)) / (max(validNumTrials) - min(validNumTrials));
        % Scale to marker size range
        markerSizes = markerSizes * (maxMarkerSize - minMarkerSize) + minMarkerSize;
    end

    % Plot percentage of reporting 'Comparison stimulus is more right' as a function of discrepancy
    scatter(validBinCenters, validPercentageRight, markerSizes, 'filled', 'DisplayName', sprintf('Standard Location %d', k));
    % Optionally, plot a line connecting the points
    plot(validBinCenters, validPercentageRight, 'LineWidth', 1, 'DisplayName', sprintf('Standard Location %d', k));

    %% sort out data needed for fitting

    D(k).standard_loc = standard_loc(k);
    D(k).comparison_loc = binCenters;
    D(k).n_r_response = numRightResponsesPerBin;
    D(k).n_trial = numTrialsPerBin;
    D(k).n_l_response = numTrialsPerBin - numRightResponsesPerBin;

end

xlim([-45 45]);
xlabel('Discrepancy relative to the standard stimulus (cm)');
ylabel('Percentage of "comparison stimulus was on the right" responses (%)');
title([ses '-2IFC discrimination']);
legend('show');
grid on;
hold off;

saveas(gcf, fullfile(out_dir, sprintf('pmf_%s_%s', str, ses)),'png');

%% fit a simple PMF
for k = 1:n_standard_locations
    % Extract data for the current standard location
    currentData = D(k);
    
    % Define the negative log-likelihood function for the current data
    nLogL = @(p) nll_gauss(p(1), p(2), p(3), currentData);
    
    % Set boundaries and initial values for the parameters
    lb = [0, -10, 0.01];
    ub = [0.06, 10, 6];
    initialization = [0, 0, 1];
    options = optimoptions(@fmincon, 'MaxIterations', 1e5, 'Display', 'off');
    
    % Fit the model using fmincon
    [est, nll] = fmincon(nLogL, initialization, [], [], [], [], lb, ub, [], options);
    
    % Display the estimated values of the parameters
    fprintf('Estimated values of the parameters for Standard Location %d:\n', k);
    disp(est);
end

%% pmf functions
function nLL = nll_gauss(lapse, mu, sigma, D)
       
    pc  = normcdf(D.comparison_loc, D.standard_loc + mu, sigma)*(1-lapse)+lapse/2;
    nLL =  -sum(log(pc.^D.n_r_response .* (1-pc).^D.n_l_response));

 end   
    
    