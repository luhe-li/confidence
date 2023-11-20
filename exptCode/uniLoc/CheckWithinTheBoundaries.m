%This function checks whether all the dot clouds are within the screen.
%To do that, we count how many clouds are off the boundaries. If NUM is
%zero, then the visual stimulus can pass this test and be presented.
function check = CheckWithinTheBoundaries(M,boxSize,ScreenInfo)
    NUM = numel(M(M<ceil(boxSize/2))) + ...
          numel(M(M([1,3],:)>(ScreenInfo.xaxis - ceil(boxSize/2)))) + ... %check x coordinates
          numel(M(M([2,4],:)>(ScreenInfo.yaxis - ceil(boxSize/2)))); %check y coordinates
    if NUM == 0; check = 1;
    else; check = 0;end
end