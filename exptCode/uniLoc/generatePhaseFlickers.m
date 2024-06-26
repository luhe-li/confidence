%This function prepares the visual stimuli by putting all the blobs together
%and convert it to a texture. Using function Screen('DrawTexture') is faster.
function flickerTx = generatePhaseFlickers(VSinfo,ExpInfo,ScreenInfo,windowPtr,i)
    %size of each check pixel is 4 by 4
    replictation = VSinfo.FlickerNumPixel; 
    height = 380;
    
    %first compute the location of the visual stimulus in pixels
    loc_pixel = round(ExpInfo.randVisPixel(i));
    targetLoc = [ScreenInfo.yaxis-ScreenInfo.liftingYaxis,ScreenInfo.xmid + loc_pixel];
    num_dots = 1000;
    RNcoordinates = round(rand(num_dots,2) .* repmat([ScreenInfo.yaxis/replictation-1,ScreenInfo.xaxis/replictation-1],num_dots,1)) + 1; % the -1 and +1 was in order to prevent index = 0;
    
    coherPDF = mvnpdf(RNcoordinates,targetLoc./replictation,[40,0;0,80]);
    coherPDF = coherPDF ./ sum(coherPDF) .*60;
    num_frames = VSinfo.numFrames;
    frameBin = randi([0 1],num_dots,num_frames);
    wavelength = 60;
    sinewave = (1-sin(linspace(0,pi,wavelength)))';
    ideal = repmat(sinewave,num_frames./wavelength,1);
    subideal = repmat(ideal,2,1);
    randInx = rand(num_dots,1)<coherPDF;
    frameBin(randInx,:) = ideal;
    outOfPhaseInd = randi([1,wavelength],sum(~randInx),1);
    frameBin(~randInx,:) = subideal(outOfPhaseInd:outOfPhaseInd+wavelength);
    flcM_repeated = NaN(ScreenInfo.yaxis,ScreenInfo.xaxis,VSinfo.GWNnumFrames);
    for n = 1:VSinfo.numFrames
        flkM = zeros(ScreenInfo.yaxis/replictation,ScreenInfo.xaxis/replictation);
        flkM(sub2ind(size(flkM),RNcoordinates(:,1), RNcoordinates(:,2))) = frameBin(:,n) .* 255;
        noise = kron(flkM(:,:),ones(replictation));
        % empty the upper and lower band of noisy background
        noise(1:(ScreenInfo.y2_lb-height),:) = 0;
        noise((ScreenInfo.y2_ub+height):end,:) = 0;
        % overlay with gray background
        flcM_repeated(:,:,n) = (noise);
%         gwnM_repeated(:,:,n) = kron(gwnM(:,:,n),ones(replictation));
        % cut the upper and lower region
        %turn it to texture
        flickerTx(n) =  Screen('MakeTexture', windowPtr, flcM_repeated(:,:,n));
    end
end