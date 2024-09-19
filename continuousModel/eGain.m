function [opt_est,opt_radius,opt_gain] = eGain(myPDF, maxScore, minScore, elbow, center_axis)
% myPDF    : two-dimensional, size equals [trial, posterior defined by center_axis]
% maxScore : the maximum possible score given to the participant
% minScore : the minimum possible score given to the participant
% elbow    : the length where the score no longer decreases and stays at
%            minScore
% center_axis : grid of possible estimates

myPDF = myPDF ./ sum(myPDF,2);
n_trial = size(myPDF,1);
n_est = numel(center_axis);
step = center_axis(2) - center_axis(1);
idx_elbow = round(elbow/step);
[maxGain,optRadius] = deal(NaN(n_trial, n_est));

% loop by possible estimate on the axis
for xx = 1:n_est

    % estimate and max confidence for each possible loctaion, same across trials
    idx_est = xx;
    idx_conf = min([idx_est - 1, n_est - idx_est]);

    % same error function across trials, only compute for nonzero max confidence radius
    if idx_conf < 1
        maxGain(:, xx) = 0;
        optRadius(:, xx) = 0;
    else

        % for each estimation location estX, it has its distinct max radius
        confRadius = 0 : idx_conf;

        % calculate the ratio between the current radius and
        % the radius where the score stops dropping (because it's too wide)
        % if confRadius >= elbow, this ratio is >= 1, participant gets the
        % minimum point
        % if confRadius < elbow, this ratio < 1, participant gets the score
        % proportional to their confRadius. e.g. ratio = 0.25, so 25% of the
        % score is deduced, they get 75% score
        lengthRatio = confRadius ./ idx_elbow;

        % This is the max possible penalty they can get. when ratio = 1, they
        % get the entire penalty.
        penaltyRange = maxScore - minScore;

        % rawScore is calculated by deducting the penalty from the max possible
        % points
        rawScore = maxScore - lengthRatio .* penaltyRange;

        % We don't let participants get negative scores, so we will take the
        % max between their raw score or the minimum score
        costFun = max(rawScore, minScore);

        % loop by trials
        for tt = 1:n_trial

            % Error pdf: For this specific estX(i) that we are working on, truncate
            % the part of the pdf from the minimum radius to the maximum radius
            % we do left and right because we can't assume symmetry in PDF around
            % estX

            erPDFright = myPDF(tt, (idx_est+1) : (idx_est + idx_conf));
            erPDFleft  = myPDF(tt, (idx_est-1) : -1 : (idx_est - idx_conf));

            % Error cdf, i.e. erf: Transform the PDF from left and right into CDF,
            % add them up along with the probability density at the estimated
            % location estX. The sum of this addition is our "error function".
            % We concatenate this with the probability density at estX because that
            % is the CDF when the confidence radius encloses nothing but the estX.
            % Conceptually, this error CDF is the probability that the stimulus is
            % in the range enclosed by each possible radii.
            % e.g. p(stimulus is within 20 pixels from estX) = erCDF(20)
            erCDF = [myPDF(tt,idx_est) , cumsum(erPDFright) + cumsum(erPDFleft) + myPDF(tt,idx_est)];

            % The expected gain is taking the dot product between the cost function
            % and the error CDF.
            % e.g. p(stimulus is within 20 pixels from estX) * cost(radius = 20)
            % erCDF(20) * costFun(20)
            gainFun = costFun .* erCDF;

            % We look for the index that corresponds to the maximum expected gain.
            % That index is the value of the radius we are looking for.
            [maxGain(tt, xx), optRadius(tt, xx)] = max(gainFun);

        end
    end

end

% convert optimal radius from axis unit to cm
optRadiusCM = optRadius.*step;

% find max gain across all possible estimated pixels for each trial
[opt_gain, idx_opt] = max(maxGain, [], 2);
opt_est = center_axis(idx_opt);
opt_radius = zeros(size(center_axis));
for tt = 1:numel(idx_opt)
    opt_radius(tt) = optRadiusCM(tt,idx_opt(tt));
end

demo_model([],[],[],[],[],[],confRadius.*step, costFun, erCDF, gainFun)

end

%% check plot

% figure; hold on
% for tt = 1:100
%     plot(myPDF(tt,:))
% end

% figure; hold on
% plot(confRadius.*step, costFun);
% plot(confRadius.*step, erCDF);
% plot(confRadius.*step, gainFun);
% xlabel('Confidence radius (cm)')

%% demo figure

% figure(f1)
% plot(confRadius.*step, costFun);
% plot(confRadius.*step, erCDF);
% plot(confRadius.*step, gainFun);
% xlabel('Confidence radius (cm)')

