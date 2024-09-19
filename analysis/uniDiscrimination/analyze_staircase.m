clear; clc; %close all;

%% set up

sub = 'LL';
ses = 'V2';
save_fig = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir, 'data','uniDiscrimination_dot')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load data
load(sprintf('uniDisDot_sub-%s_ses-%s', sub, ses));
n_staircase = ExpInfo.n_staircase;

% replace nan trials
Resp.discrepancy(1,39:40) = min(ExpInfo.comparison_loc);
Resp.discrepancy(2,39:40) = max(ExpInfo.comparison_loc);
Resp.correct(1,39:40) = 1;
Resp.correct(2,39:40) = 1;


%% Plot the interleaved staircases

% Prepare figure
figure;
hold on;

trialNum = 1:ExpInfo.n_trial;
n_easy_trials = ExpInfo.n_easy_trial_per_s;
speaker_cm = ExpInfo.speaker_level_cm;

% For Condition 1
condition1 = 1; % Assuming condition 1 is in row 1
comparison_loc1 = Resp.comparison_loc(condition1, 1:(end - n_easy_trials));
correct1 = Resp.correct(condition1, 1:(end - n_easy_trials));

% Determine comparison values based on 'ses'
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

% For Condition 2
condition2 = 2; % Assuming condition 2 is in row 2
comparison_loc2 = Resp.comparison_loc(condition2, 1:(end - n_easy_trials));
correct2 = Resp.correct(condition2, 1:(end - n_easy_trials));

% Determine comparison values based on 'ses'
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
ylim([-80 80])
title([ses '-2IFC Discrimination']);

legend([h_line1, h_line2], '1-up-2-down', '2-up-1-down');
grid on;
hold off;


%% plot raw data for PMF
% Flatten the matrices into vectors and remove NaN values
Discrepancies = Resp.discrepancy(:);
Responses = Resp.resp(:);
validIndices = ~isnan(Discrepancies) & ~isnan(Responses);
Discrepancies = Discrepancies(validIndices);
Responses = Responses(validIndices);

% Determine bin centers and edges
if strcmp(ses, 'A')
    % Convert speaker index to cm
    comparison_loc = speaker_cm(ExpInfo.comparison_loc);
    Discrepancies = Discrepancies * (speaker_cm(2) - speaker_cm(1));
    bin_width = 5;
else
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

% Bin the discrepancies
[~, ~, binIdx] = histcounts(Discrepancies, binEdges);

% Initialize arrays for counts
numTrialsPerBin = zeros(size(binCenters));
numRightResponsesPerBin = zeros(size(binCenters));

% Loop over bins to count responses
for b = 1:length(binCenters)
    idx = binIdx == b;
    numTrialsPerBin(b) = sum(idx);
    numRightResponsesPerBin(b) = sum(Responses(idx) == 1);
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
figure;
scatter(validBinCenters, validPercentageRight, markerSizes, 'ko', 'filled');
hold on;
% Optionally, plot a line connecting the points
plot(validBinCenters, validPercentageRight, 'k-', 'LineWidth', 1);

xlim([-45 45]);
xlabel('Discrepancy relative to the standard stimulus (cm)');
ylabel('Percentage of "comparison stimulus was on the right" responses (%)');
title([ses '-2IFC discrimination']);
grid on;
hold off;

