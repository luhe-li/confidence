function updated_loc = findNextLocation(i,j,current_loc, current_resp, history, ExpInfo)
    
    condition = ExpInfo.condition(j, i);
    TrialNumber = i; % Calculate trial number based on history

    % Condition 1: 1-up-2-down, target starts on the LEFT side
    if rem(condition, 2) == 1  
        if TrialNumber == 1  % If this is the first trial
            if current_resp == -1  % Response is "LEFT", which is correct
                updated_loc = current_loc + ExpInfo.StepSizes(1);  % Make it harder
            else
                updated_loc = current_loc;  % Keep difficulty the same
            end
        else
            % Find consecutive "RIGHT" responses in history + current_resp
            numPairOfRight = findPair([history, current_resp], [1 1 -1]);

            % Adjust step size based on the number of consecutive "RIGHT" responses
            if numPairOfRight == 0
                stepsize = ExpInfo.StepSizes(1);
            elseif numPairOfRight == 1
                stepsize = ExpInfo.StepSizes(2);
            else
                stepsize = ExpInfo.StepSizes(3);
            end

            if current_resp == -1 && history(end) == -1
                updated_loc = current_loc + stepsize;  % Make it harder
                disp('Correct twice, decrease distance')
            elseif current_resp == 1
                updated_loc2 = current_loc - stepsize;  % Make it easier
                disp('Wrong once, widen distance')
            else
                updated_loc = current_loc;  % Keep difficulty the same
                disp('Correct once, keep distance')
            end
        end
    else  % Condition 2: 1-down-2-up, target starts on the RIGHT side
        if TrialNumber == 1  % If this is the first trial
            if current_resp == 1  % Response is "RIGHT", which is correct
                updated_loc = current_loc - ExpInfo.StepSizes(1);  % Make it easier
            else
                updated_loc = current_loc;  % Keep difficulty the same
            end
        else
            % Find consecutive "LEFT" responses in history + current_resp
            numPairOfLeft = findPair([history current_resp], [-1 -1 1]);

            % Adjust step size based on the number of consecutive "LEFT" responses
            if numPairOfLeft == 0
                stepsize = ExpInfo.StepSizes(1);
            elseif numPairOfLeft == 1
                stepsize = ExpInfo.StepSizes(2);
            else
                stepsize = ExpInfo.StepSizes(3);
            end

            if current_resp == -1 && history(end) == -1
                updated_loc = current_loc + stepsize;  % Make it easier
                disp('Wrong twice, widen distance')
            elseif current_resp == 1
                updated_loc = current_loc - stepsize;  % Make it harder
                disp('Correct, decrease distance')
            else
                updated_loc = current_loc;  % Keep difficulty the same
                disp('Wrong once, keep distance')
            end

        end
    end
end

function numPair = findPair(history, number)
    length_h = length(history);
    length_n = length(number);
    test_trials = length_h - (length_n - 1);
    numPair = 0;
    for i = 1:test_trials
        if isequal(history(i:(i + length_n - 1)), number)
            numPair = numPair + 1;
        end
    end
end