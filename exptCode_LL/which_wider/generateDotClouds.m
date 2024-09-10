%This function prepares the visual stimuli by putting all the blobs together
%and convert it to a texture. Using function Screen('DrawTexture') is faster.
function dotClouds = generateDotClouds(windowPtr,locations,VSinfo,ScreenInfo)
    locations = round(locations);
    for i = 1:VSinfo.n_dots 
        %For each centroid, we put a square centered around it. Within this
        %box, we draw a gaussian distribution and add it to a grey
        %background. Keep doing this for each centroid.
        VSinfo.blankScreen((locations(1,i)-floor(VSinfo.boxSize/2)):...
            (locations(1,i)+floor(VSinfo.boxSize/2)),(locations(2,i)-...
            floor(VSinfo.boxSize/2)):(locations(2,i)+...
            floor(VSinfo.boxSize/2))) = VSinfo.Cloud;
        VSinfo.greyScreen = VSinfo.greyScreen + VSinfo.blankScreen; 
        %keep adding the cloud to the grey background
        VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis); 
        %restore the blank screen
    end     
    %if the color value hits beyond 255 (maximum), then we clip it
    VSinfo.greyScreen(VSinfo.greyScreen>255)=255;
    %Turn the matrix to texture    
    dotClouds =  Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
end
