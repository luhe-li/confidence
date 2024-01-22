%This function generates visual masker, which is 900 ms long (54 frames).
function gwnBackground = generateNoisyBackground(VSinfo,ScreenInfo,windowPtr)
    %size of each check pixel is 4 by 4
    replictaion = VSinfo.GWNnumPixel; 
    %randomly generate brightness between 0 and 255. However, we don't
    %generae 54 frames. We only generate 10 frames and use them repeatedly
    gwnM        = 255.*rand(ScreenInfo.yaxis/replictaion,ScreenInfo.xaxis/...
                    replictaion,VSinfo.GWNnumFrames);
    %initialize
    gwnM_repeated = NaN(ScreenInfo.yaxis,ScreenInfo.xaxis,VSinfo.GWNnumFrames);
    for n = 1:VSinfo.GWNnumFrames
        gwnM_repeated(:,:,n) = kron(gwnM(:,:,n),ones(replictaion));
        %turn it to texture
        gwnBackground(n) =  Screen('MakeTexture', windowPtr, gwnM_repeated(:,:,n));
    end
end