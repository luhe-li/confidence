function [opt_radius,opt_gain, func] = eGain_MAP(myPDF, estX, maxScore, minScore, elbow, center_axis)
% myPDF    : two-dimensional, size equals [trial, posterior defined by center_axis]
% estX     : the estimated location for each trial, size equals [trial, 1]
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
[opt_gain,idx_opt_radius] = deal(NaN(n_trial,1));
func = repmat(struct('costFun', [], 'erCDF', [], 'gainFun', []), n_trial, 1);

% Calculate the closest elements for each value in estX and indices
[~, vec_idx_est] = min(abs(center_axis - estX), [], 2);

% picking the minimum allowed value as the max possible radius length
% estX - screen_cm/2 is how far the estimation point is from the left edge of screen
% screen_cm/2 - estX is how far from the right edge of screen
% if confRadius is more than the min of these two, it would grow out of the
% screen
vec_idx_maxconf = min([vec_idx_est - 1, n_est - vec_idx_est],[],2);

% loop by trials
for tt = 1:length(estX)

    % index of approxiamted estimate on the sampling axis of this trial
    idx_est = vec_idx_est(tt);
    idx_conf = vec_idx_maxconf(tt);

    if idx_conf < 1
        opt_gain(tt) = 0;
        idx_opt_radius(tt) = 0;
        func(tt).costFun = [];
        func(tt).erCDF = [];
        func(tt).gainFun = [];
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
        func(tt).costFun = costFun;
        func(tt).erCDF = erCDF;
        func(tt).gainFun = gainFun;

        % We look for the index that corresponds to the maximum expected gain.
        % That index is the radius we are looking for in the unit of axis defined.
        [opt_gain(tt), idx_opt_radius(tt)] = max(gainFun);

    end
end

% convert optimal radius from axis unit to cm
opt_radius = idx_opt_radius.*step;

% ii=18;
% temp = 0:vec_idx_maxconf(ii);
% demo_model2(myPDF(ii,:), center_axis, estX(ii), opt_radius(ii), temp.*step,...
%     func(ii).costFun, func(ii).erCDF, func(ii).gainFun)

end

% % check posterior
% 
% figure; hold on
% for tt = 1:2
%     plot(myPDF(tt,:))
% end
% 
% % check functions
% figure; hold on
% plot(confRadius.*step, costFun);
% plot(confRadius.*step, erCDF);
% plot(confRadius.*step, gainFun);
% xlabel('Confidence radius (cm)')