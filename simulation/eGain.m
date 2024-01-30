function [optRadius,expGain] = eGain(myPDF, estX, maxPoint, minPoint, elbow, screenX)
% myPDF : one-dimensional, size equals [1,screenX]
% estX  : the estimated location
% elbow : the point where the cost hits minPoint
% screenX : ScreenInfo.xaxis

myPDF = myPDF ./ sum(myPDF);

confRangeMax = min( [estX - 1, screenX - estX] );
confRange = 0 : confRangeMax;
costFun = max( maxPoint - confRange ./ elbow .* (maxPoint - minPoint) , minPoint);

erPDF = myPDF((estX+1) : (estX + confRangeMax)); % error pdf
erCDF = [myPDF(estX) , cumsum(erPDF) .* 2 + myPDF(estX)]; % error cdf, i.e. erf

expGain = costFun .* erCDF;
[~,optRadius] = max(expGain);
end