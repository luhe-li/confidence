%Since the visual stimulus is generated randomly given a predetermined
%test stimulus location, the actual center of the 10 Gaussian blobs will 
%be different from the predetermined test location. This function moves the 
%whole collection of 10 blobs to align with the predetermined test location.
function locations = shiftDotClouds(currentLocations,centroid_predetermined,...
    ScreenInfo)
    %first initialize
    locations      = NaN(size(currentLocations));
    %get the actual center of the 10 Gaussian blobs
    averagedXY     = mean(currentLocations,2);
    centroid_xAxis = averagedXY(1);
    %compute the difference
    diff           = ScreenInfo.xmid + ...
                        centroid_predetermined - centroid_xAxis;
    %shifting the horizontal location of the blobs to match the
    %predtermined test location
    locations(1,:) = currentLocations(1,:) + diff;
    %we don't adjust the vertical location of the stimulus
    locations(2,:) = currentLocations(2,:);
end