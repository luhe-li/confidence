clear; close all;

sub_init = 'OY';
session = '1';
exp = 'biLoc';

% manage path
[project_dir, ~] = fileparts(pwd);
if ~strcmp(fileparts(project_dir), 'confidence')
    [project_dir, ~] = fileparts(project_dir);
end
data_dir     = fullfile(project_dir, 'data', exp);

file_path = fullfile(data_dir, sprintf('biLoc_sub-%s_ses-%s_incomplete.mat', sub_init, session));
load(file_path);


% Check if ExpInfo.nTrials matches the length of Resp
if ExpInfo.nTrials ~= length(Resp)
    % Calculate the number of complete sets of trials
    nCompleteSets = floor(length(Resp) / ExpInfo.nLevel);

    % Adjust Resp to only include complete sets of trials
    Resp = Resp(1:nCompleteSets * ExpInfo.nLevel);

    % Update ExpInfo
    ExpInfo.nRep = nCompleteSets;
    ExpInfo.nTrials = length(Resp);

    % Update variables including 'rand' in their names
    randFields = fieldnames(ExpInfo);
    for i = 1:length(randFields)
        if contains(randFields{i}, 'rand')
            ExpInfo.(randFields{i}) = ExpInfo.(randFields{i})(1:ExpInfo.nTrials);
        end
    end
end

% Sort trials by post-cue first (A=1, V=2), and location level next
sortMatrix = [[Resp.post_cue]; [Resp.target_idx]]';
[~, order] = sortrows(sortMatrix);
sortedResp = Resp(order);

% Save the updated data
outFileName = sprintf('biLoc_sub-%s_ses-%s.mat', sub_init, session);
save(fullfile(data_dir, outFileName), 'Resp', 'sortedResp', 'ExpInfo', 'ScreenInfo', 'VSinfo', 'AudInfo');
