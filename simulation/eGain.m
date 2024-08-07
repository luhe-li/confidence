function [optRadius,expGain] = eGain(myPDF, estX, maxScore, minScore, elbow, screenX)
% myPDF    : one-dimensional, size equals [1,screenX]
% estX     : the estimated location
% maxScore : the maximum possible score given to the participant
% minScore : the minimum possible score given to the participant
% elbow    : the length where the score no longer decreases and stays at
%            minScore
% screenX  : ScreenInfo.xaxis

myPDF = myPDF ./ sum(myPDF,2); % normalize the PDF before any operation

confRadiusMax = min( [estX - 1, screenX - estX],[],2); 
% picking the minimum allowed value as the max possible radius length
% estX - 1 is how far the estimation point is from the left edge of screen
% screenX - estX is how far from the right edge of screen
% if confRadius is more than the min of these two, it would grow out of the
% screen

optRadius = NaN(length(estX),1);
% initialize an array to store the best radius for each estimation location
% estX

for i = 1:length(estX)
    
    confRadius = 0 : confRadiusMax(i); 
    % for each estimation location estX, it has its distinct max radius

    lengthRatio = confRadius ./ elbow;
    % calculate the ratio between the current radius and 
    % the radius where the score stops dropping (because it's too wide) 
    % if confRadius >= elbow, this ratio is >= 1, participant gets the
    % minimum point
    % if confRadius < elbow, this ratio < 1, participant gets the score
    % proportional to their confRadius. e.g. ratio = 0.25, so 25% of the 
    % score is deduced, they get 75% score
    
    penaltyRange = maxScore - minScore;
    % This is the max possible penalty they can get. when ratio = 1, they
    % get the entire penalty. 
    rawScore = maxScore - lengthRatio .* penaltyRange;
    % rawScore is calculated by deducting the penalty from the max possible
    % points
    costFun = max(rawScore, minScore);
    % We don't let participants get negative scores, so we will take the
    % max between their raw score or the minimum score

    erPDFright = myPDF(i,(estX(i)+1) : (estX(i) + confRadiusMax(i))); 
    erPDFleft  = myPDF(i,(estX(i)-1) : -1 : (estX(i) - confRadiusMax(i))); 
    % Error pdf: For this specific estX(i) that we are working on, truncate
    % the part of the pdf from the minimum radius to the maximum radius
    % we do left and right because we can't assume symmetry in PDF around
    % estX
    erCDF = [myPDF(i,estX(i)) , cumsum(erPDFright) + cumsum(erPDFleft) + myPDF(i,estX(i))]; 
    % Error cdf, i.e. erf: Transform the PDF from left and right into CDF,
    % add them up along with the probability density at the estimated
    % location estX. The sum of this addition is our "error function". 
    % We concatenate this with the probability density at estX because that
    % is the CDF when the confidence radius encloses nothing but the estX. 
    % Conceptually, this error CDF is the probability that the stimulus is
    % in the range enclosed by each possible radii. 
    % e.g. p(stimulus is within 20 pixels from estX) = erCDF(20)
    eGain = costFun .* erCDF;
    % The expected gain is taking the dot product between the cost function
    % and the error CDF. 
    % e.g. p(stimulus is within 20 pixels from estX) * cost(radius = 20)
    % = erCDF(20) * costFun(20)
    [~,optRadius(i)] = max(eGain);
    % We look for the index that corresponds to the maximum expected gain. 
    % That index is the value of the radius we are looking for.
    expGain(i).eGain = eGain;
    % We also save this entire eGain array in case we need it later
    % (totally fine to remove, just remember to change function output)
end

end





