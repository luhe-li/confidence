function leaderboardText = updateLeaderboardBimodal(outDir, ExpInfo, Resp)
    % updateLeaderboardUnimodal Updates or creates the leaderboard for a unimodal psychology task.
    %
    % Syntax:
    %   leaderboardText = updateLeaderboardUnimodal(outDir, ExpInfo, Resp)
    %
    % Inputs:
    %   outDir    - String. The output directory where 'sub_info.csv' is stored.
    %   ExpInfo   - Struct with fields:
    %                .subjInit - String. Subject identifier.
    %                .session  - Char. Session type ('A' or 'V').
    %   Resp      - Struct with field:
    %                .point    - Numeric array. Points earned in the session.
    %
    % Outputs:
    %   leaderboardText - String. Formatted leaderboard text for the current session.

    % Define the CSV file path
    csvPath = fullfile(outDir, 'sub_info.csv');

    % Check if sub_info.csv exists
    if ~exist(csvPath, 'file')
        % Create a new table with appropriate variable types
        sub_info = table('Size', [0 4], ...
                         'VariableTypes', {'cellstr', 'double', 'double', 'double'}, ...
                         'VariableNames', {'sub_name', 'sub_id', 'ses_a_pt', 'ses_v_pt'});
    else
        % Read the existing CSV file
        sub_info = readtable(csvPath);
    end

    % Check if the current subject exists in the table
    subjectIdx = find(strcmp(sub_info.sub_name, ExpInfo.subjInit), 1);

    if isempty(subjectIdx)
        % Subject does not exist; create a new entry

        % Assign a new sub_id
        if isempty(sub_info.sub_id)
            new_sub_id = 1;
        else
            new_sub_id = max(sub_info.sub_id) + 1;
        end

        % Calculate total points earned in the session
        total_pt = sum([Resp.point]);

        % Initialize session points
        new_ses_a_pt = 0;
        new_ses_v_pt = 0;

        % Assign points based on the session type
        if strcmpi(ExpInfo.session, 'A')
            new_ses_a_pt = total_pt;
        elseif strcmpi(ExpInfo.session, 'V')
            new_ses_v_pt = total_pt;
        else
            error('Invalid session type. Must be ''A'' or ''V''.');
        end

        % Create a new row for the subject
        newRow = table({ExpInfo.subjInit}, new_sub_id, new_ses_a_pt, new_ses_v_pt, ...
            'VariableNames', {'sub_name', 'sub_id', 'ses_a_pt', 'ses_v_pt'});

        % Append the new row to the existing table
        sub_info = [sub_info; newRow];

        % Write the updated table back to the CSV file
        writetable(sub_info, csvPath);

        % Display confirmation message
        disp('New subject added and points recorded for this session.');
    else
        % Subject exists; update the existing entry

        % Calculate total points earned in the session
        total_pt = sum([Resp.point]);

        % Update points based on the session type
        if strcmpi(ExpInfo.session, 'A')
            sub_info.ses_a_pt(subjectIdx) = sub_info.ses_a_pt(subjectIdx) + total_pt;
        elseif strcmpi(ExpInfo.session, 'V')
            sub_info.ses_v_pt(subjectIdx) = sub_info.ses_v_pt(subjectIdx) + total_pt;
        else
            error('Invalid session type. Must be ''A'' or ''V''.');
        end

        % Write the updated table back to the CSV file
        writetable(sub_info, csvPath);

        % Display confirmation message
        disp('Existing subject updated with new points for this session.');
    end

    % Initialize the output text
    leaderboardText = '';

    % Generate leaderboard text based on the current session
    if strcmpi(ExpInfo.session, 'A')
        % Sort the table based on auditory points in descending order
        sorted_table = sortrows(sub_info, 'ses_a_pt', 'descend');

        % Assign ranks based on sorted order
        sorted_table.rank = (1:height(sorted_table))';

        % Initialize the leaderboard string with an explanation
        explanation = 'Leaderboard for Unimodal Auditory Localization Task:\n\n';

        % Create a header for the table
        header = sprintf('%-20s %-5s %-10s', 'Name', 'Rank', 'Points');

        % Start building the leaderboard text
        leaderboardText = [explanation, header, '\n'];

        % Iterate through each row to append subject details
        for i = 1:height(sorted_table)
            name = sorted_table.sub_name{i};
            rank = sorted_table.rank(i);
            points = sorted_table.ses_a_pt(i);
            rowText = sprintf('%-20s %-5d %.2f', name, rank, points);
            leaderboardText = [leaderboardText, rowText, '\n'];
        end

    elseif strcmpi(ExpInfo.session, 'V')
        % Sort the table based on visual points in descending order
        sorted_table = sortrows(sub_info, 'ses_v_pt', 'descend');

        % Assign ranks based on sorted order
        sorted_table.rank = (1:height(sorted_table))';

        % Initialize the leaderboard string with an explanation
        explanation = 'Leaderboard for Audiovisual Localization Task:\n\n';

        % Create a header for the table
        header = sprintf('%-20s %-5s %-10s', 'Name', 'Rank', 'Points');

        % Start building the leaderboard text
        leaderboardText = [explanation, header, '\n'];

        % Iterate through each row to append subject details
        for i = 1:height(sorted_table)
            name = sorted_table.sub_name{i};
            rank = sorted_table.rank(i);
            points = sorted_table.ses_v_pt(i); % Corrected to ses_v_pt
            rowText = sprintf('%-20s %-5d %.2f', name, rank, points);
            leaderboardText = [leaderboardText, rowText, '\n'];
        end

    else
        % For any other session types, return an empty string
        leaderboardText = '';
    end
end