function [H, P, STATS,d] = af_cohenD(X1,X2)
% calculate cohens d for ttest and such
% d = M1-M2/SDpooled
% SD pooled = sqrt(SD1^2+SD2^2/2)
if nargin < 2 %one vector is being compared to 0 see if it's significantly different
    % calculate mean
    meanX1 = mean(X1);
    % calculated std
    stdX1 = std(X1);
    % calculate Cohen's d
    d = (meanX1-0)/stdX1;
    % also do a ttest bcuz why not - gonna need it anyway
    [H,P,~,STATS] = ttest(X1,0,'Tail','both','Alpha',0.05);
elseif nargin == 2 %two vectors are being compared to see if they are significantly different from each other
    % calculate mean
    meanX1 = mean(X1);
    meanX2 = mean(X2);
    % calculated std
    stdX1 = std(X1);
    stdX2 = std(X2);
    % calculated the pooled std
    SD_pooled = sqrt((stdX1^2+stdX2^2)/2);
    % calculate Cohen's d
    d = (meanX1-meanX2)/SD_pooled;
    % also do a ttest bcuz why not - gonna need it anyway
    [H,P,~,STATS] = ttest2(X1,X2,'Vartype','equal','Tail','both'); 
else
    disp('Error!');
end

