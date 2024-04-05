function [data_by_diff, raw_diffs] = org_by_raw_diffs_4D(org_data, sA)

sV = sA;
raw_diffs = unique(sA' - sV)';
num_diffs = length(raw_diffs);
data_by_diff = nan(num_diffs, 2, 2); % diff x cue x reliability

for i = 1:num_diffs
    
    diff = raw_diffs(i);
    
    % Find pairs of audIdx and visIdx that match this difference
    [audPairs, visPairs] = find((sA' - sV) == diff);
    
    tempData = [];
    % For each pair, extract and store the corresponding data
    for j = 1:numel(audPairs)
        
        % Extract data for this specific audIdx and visIdx pair across all other dimensions
        tempData = cat(3, tempData, squeeze(org_data(audPairs(j), visPairs(j), :, :)));
    
    end
    
    % take mean over the repeats if conditions that satisfied the discrepancy > 1
    if j > 1
        m_temp_data = mean(tempData,3);
    else
        m_temp_data = tempData;
    end
    
    % store organized data by discrepancy
    data_by_diff(i,:,:) = m_temp_data;
        
end

end