% This function organizes 4D auditory and visual data based on their *raw* discrepancies.
% It computes the differences between every pair of auditory and visual indices,
% groups the data by these differences, and *averages* the data for each group
% if there are multiple entries per group.

%-------------------------------------------------------------------------------
% Inputs:
% org_data: The original 4D dataset to be organized, with dimensions representing
%           [auditory location, visual location, cue, reliability.
% sA: A vector representing the auditory location indices. The subtraction of
%     sA from itself is used to find the raw differences.

% Outputs:
% data_by_diff: A 3D array organized by discrepancy, cue, and reliability dimensions.
%               The first dimension represents unique discrepancies calculated from sA.
%               The size is [number of unique discrepancies, 2, 2].
% raw_diffs: A vector containing the unique raw discrepancies calculated from sA.
%-------------------------------------------------------------------------------

function [data_by_diff, raw_diffs] = org_by_raw_diffs_4D(data, sA)

sV = sA;
raw_diffs = unique(round((sV' - sA),2))';
num_diffs = length(raw_diffs);
data_by_diff = nan(num_diffs, 2, 2); % diff x cue x reliability

for i = 1:num_diffs
    
    diff = raw_diffs(i);
    
    % Find pairs of audIdx and visIdx that match this difference
    [visPairs, audPairs] = find(round((sV' - sA),2) == diff);
    
    tempData = [];
    % For each pair, extract and store the corresponding data
    for j = 1:numel(audPairs)
        
        % Extract data for this specific audIdx and visIdx pair across all other dimensions
        tempData = cat(3, tempData, squeeze(data(audPairs(j), visPairs(j), :, :)));
    
    end
    
    % Average data across locations for the current discrepancy, if applicable
    if j > 1
        m_temp_data = mean(tempData,3);
    else
        m_temp_data = tempData;
    end
    
    % Store the organized (and possibly averaged) data by discrepancy
    data_by_diff(i,:,:) = m_temp_data;
        
end

end