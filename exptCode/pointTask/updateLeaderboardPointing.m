function leaderboardText = updateLeaderboardPointing(outDir, ExpInfo, Resp)
    % updateLeaderboard Updates or creates the leaderboard for a psychology task.
    %
    % Syntax:
    %   leaderboardText = updateLeaderboard(projectDir, ExpInfo, Resp)
    %
    % Inputs:
    %   projectDir - String. The root directory of the project.
    %   ExpInfo    - Struct with fields:
    %                .subjInit - String. Subject identifier.
    %                .session  - Char. Session type ('A' or 'V').
    %   Resp       - Struct with field:
    %                .point    - Numeric array. Points earned in the session.
    %
    % Outputs:
    %   leaderboardText - String. Formatted leaderboard text for session 'A'.
    %                     Empty string for other sessions.
    
    % Define the CSV file path
    csvPath = fullfile(outDir, 'sub_info.csv');
    
    % Check if sub_info.csv exists
    if ~exist(csvPath, 'file')
        % Create a new table with appropriate variable types
        sub_info = table('Size', [0 3], ...
                         'VariableTypes', {'cellstr', 'double', 'double'}, ...
                         'VariableNames', {'sub_name', 'sub_id', 'pt'});
    else
        % Read the existing CSV file
        sub_info = readtable(csvPath);
    end
    
    % Check if the current subject exists in the table
    subjectExists = any(strcmp(sub_info.sub_name, ExpInfo.subjInit));
    
    if ~subjectExists
        % Assign a new sub_id
        if isempty(sub_info.sub_id)
            new_sub_id = 1;
        else
            new_sub_id = max(sub_info.sub_id) + 1;
        end
        
        % Calculate total points earned in the session
        new_pt = sum([Resp.point]);

        % Create a new row for the subject
        newRow = table({ExpInfo.subjInit}, new_sub_id, new_pt, ...
            'VariableNames', {'sub_name', 'sub_id', 'pt'});
        
        % Append the new row to the existing table
        sub_info = [sub_info; newRow];
        
        % Write the updated table back to the CSV file
        writetable(sub_info, csvPath);
        
        % Display confirmation message
        disp('Point saved to record for this session.');
    else
        % Subject already exists; you can implement updating logic here if needed
        disp('Subject already exists. No new entry added.');
    end
    
    % Initialize the output text
    leaderboardText = '';
    
    % Read the updated sub_info table
    sub_info = readtable(csvPath);
    
    % Sort the table based on auditory points in descending order
    sorted_table = sortrows(sub_info, 'pt', 'descend');
    
    % Assign ranks based on sorted order
    sorted_table.rank = (1:height(sorted_table))';
    
    % Initialize the leaderboard string with an explanation
    explanation = 'Leaderboard for the Pointing Task:\n\n';
    
    % Create a header for the table
    header = sprintf('%-20s %-5s %-10s', 'Name', 'Rank', 'Points');
    
    % Start building the leaderboard text
    leaderboardText = [explanation, header, '\n'];
    
    % Iterate through each row to append subject details
    for i = 1:height(sorted_table)
        name = sorted_table.sub_name{i};
        rank = sorted_table.rank(i);
        points = sorted_table.pt(i);
        rowText = sprintf('%-20s %-5d %.2f', name, rank, points);
        leaderboardText = [leaderboardText, rowText, '\n'];
    end
end