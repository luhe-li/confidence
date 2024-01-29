
for sub = 1:5
    flnm = sprintf('df_overall_S%i.csv',sub);
    df = readtable(flnm);
    data = table2struct(df);

    a = [data.StimulusLocationA; data.ResponsesA]';
    v = [data.StimulusLocationV; data.ResponsesV; data.sesN]';


    %%
    % Sort 'a' by 'StimulusLocationA' (first column)
    [~, idx] = sort(a(:, 1));
    sortedA = a(idx, :);

    % Find unique values of 'StimulusLocationA' and their indices
    [uniqueLocations, ~, groupIndices] = unique(sortedA(:, 1));

    % Initialize an array to store the standard deviations
    a_stdDeviations = zeros(size(uniqueLocations));

    % Calculate standard deviation for each group
    for i = 1:length(uniqueLocations)
        a_stdDeviations(i) = std(sortedA(groupIndices == i, 2));
    end

    % Display the standard deviations
    disp(a_stdDeviations);
    disp(mean(a_stdDeviations));


    %%
    % Sort 'v' by 'sesN' (third column)
    [~, idx] = sort(v(:, 3));
    v_sorted_sesN = v(idx, :);

    % Find unique values of 'sesN' and their indices
    uniqueSesN = unique(v_sorted_sesN(:, 3));

    % Initialize a matrix to store the standard deviations
    stdDeviations = [];

    % Iterate over each unique 'sesN'v
    for i = 1:length(uniqueSesN)
        % Extract rows corresponding to the current 'sesN'
        rows = v_sorted_sesN(v_sorted_sesN(:, 3) == uniqueSesN(i), :);

        % Sort these rows by 'StimulusLocationA'
        [~, idx] = sort(rows(:, 1));
        sortedRows = rows(idx, :);

        % Find unique values of 'StimulusLocationA' in this subset
        uniqueLocations = unique(sortedRows(:, 1));

        % Calculate standard deviation for each 'StimulusLocationA' within this 'sesN'
        for j = 1:length(uniqueLocations)
            stdDev = std(sortedRows(sortedRows(:, 1) == uniqueLocations(j), 2));
            stdDeviations = [stdDeviations; uniqueSesN(i), uniqueLocations(j), stdDev];
        end
    end

    % Display the standard deviations
    % disp(stdDeviations);

    %%
    % Initialize an array to store the mean standard deviations for each sesN
    meanStdDeviations = [];

    % Iterate over each unique 'sesN'
    for i = 1:length(uniqueSesN)
        % Extract standard deviations for the current 'sesN'
        sesNStdDeviations = stdDeviations(stdDeviations(:, 1) == uniqueSesN(i), 3);

        % Calculate the mean standard deviation for this 'sesN'
        meanStd = mean(sesNStdDeviations);

        % Append to the result array
        meanStdDeviations = [meanStdDeviations; uniqueSesN(i), meanStd];
    end

    % Display the mean standard deviations for each sesN
    disp(meanStdDeviations);


    %% plot
    figure
    bar(sort(meanStdDeviations(:,2)))
    l = yline(mean(a_stdDeviations),'LineWidth',2);
    xlabel('visual reliability conditions')
    ylabel('visual response s.d.')
    legend(l, 'auditory response s.d.','Location','bestoutside')
    title(sprintf('sub%i',sub))
    
    saveas(gcf, sprintf('sub%i',sub),'png')

end