function [err, y_pred] = sqErr_locallyLinearFunc(slope, intercept, x, y)
%to answer reviewer 2's question
%"I wonder whether the relation between the PSEs and the location of the 
%standard visual stimulus is indeed a linear one; for instance, looking a 
%Figure 2B (right panel) one may also consider such a relation as locally 
%linear with a slope of 1 (plus an offset) between the two central points, 
%and shallower towards the two extreme points."
    y_pred    = NaN(size(x));
    y_pred(2) = x(2) + intercept;
    y_pred(1) = y_pred(2) - slope*(x(2)-x(1));
    y_pred(3) = x(3) + intercept;
    y_pred(4) = y_pred(3) + slope*(x(4)-x(3));
    err       = sum((y - y_pred).^2);
end