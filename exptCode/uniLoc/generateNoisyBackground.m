%This function generates visual masker, which is 900 ms long (54 frames).
function gwnBackground = generateNoisyBackground(VSinfo,ScreenInfo,windowPtr)
    %size of each check pixel is 4 by 4
    replictation = VSinfo.GWNnumPixel; 
    height = 80;
    %randomly generate brightness between 0 and 255
    gwnM        = 0.3.*255.*rand(ScreenInfo.yaxis/replictation,...
        ScreenInfo.xaxis/replictation,VSinfo.GWNnumFrames);
    %initialize
    gwnM_repeated = NaN(ScreenInfo.yaxis,ScreenInfo.xaxis,VSinfo.GWNnumFrames);
    for n = 1:VSinfo.GWNnumFrames
        noise = kron(gwnM(:,:,n),ones(replictation));
        % empty the upper and lower band of noisy background
        noise(1:(ScreenInfo.y2_lb-height),:) = 0;
        noise((ScreenInfo.y2_ub+height):end,:) = 0;
        % overlay with gray background
        gwnM_repeated(:,:,n) = VSinfo.greyScreen' + noise;
%         gwnM_repeated(:,:,n) = kron(gwnM(:,:,n),ones(replictation));
        % cut the upper and lower region
        %turn it to texture
        gwnBackground(n) =  Screen('MakeTexture', windowPtr, gwnM_repeated(:,:,n));
    end
end