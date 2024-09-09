% This function sorts the 5D data by absolute discrepancies. Usually used
% when organizing bimodal confidence data.

function [conf_by_diff, all_diffs] = org_by_diffs(org_data, sA)

sV = sA;
all_diffs = unique(round(abs(sV' - sA),2))';
num_diffs = length(all_diffs);
conf_by_diff = cell(1,num_diffs);

for i = 1:num_diffs
    
    diff = all_diffs(i);
    
    % Find pairs of audIdx and visIdx that match this difference
    [audPairs, visPairs] = find(round(abs(sV' - sA),2) == diff);
    
    tempData = [];
    % For each pair, extract and store the corresponding data
    for j = 1:numel(audPairs)
        
        % Extract data for this specific audIdx and visIdx pair across all other dimensions
        tempData = cat(3, tempData, squeeze(org_data(audPairs(j), visPairs(j), :, :, :)));
    
    end
    
    % store organized data by discrepancy
    conf_by_diff{i} = tempData;
        
end

end