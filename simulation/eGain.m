function [optRadius,expGain] = eGain(myPDF, estX, maxPoint, minPoint, elbow, screenX)
% myPDF   : one-dimensional, size equals [1,screenX]
% estX    : the estimated location
% elbow   : the point where the cost hits minPoint
% screenX : ScreenInfo.xaxis

myPDF = myPDF ./ sum(myPDF,2);

confRangeMax = min( [estX - 1, screenX - estX],[],2);
optRadius = NaN(length(estX),1);
for i = 1:length(estX)
    confRange = 0 : confRangeMax(i);
    costFun = max( maxPoint - confRange ./ elbow .* (maxPoint - minPoint) , minPoint);

    erPDF = myPDF(i,(estX(i)+1) : (estX(i) + confRangeMax(i))); % error pdf
    erCDF = [myPDF(i,estX(i)) , cumsum(erPDF) .* 2 + myPDF(i,estX(i))]; % error cdf, i.e. erf

    eGain = costFun .* erCDF;
    [~,optRadius(i)] = max(eGain);
    expGain(i).eGain = eGain;
end

end