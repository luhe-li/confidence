%This function prepares the visual stimuli by putting all the blobs together
%and convert it to a texture. Using function Screen('DrawTexture') is faster.
function flickerTx = generateBinFlickers(VSinfo,ExpInfo,ScreenInfo,windowPtr,i)
    %size of each check pixel is 4 by 4
    replictation = VSinfo.FlickerNumPixel; 
    height = 380;
    
    %first compute the location of the visual stimulus in pixels
    loc_pixel = round(ExpInfo.randVisPixel(i));
    targetLoc = [ScreenInfo.yaxis-ScreenInfo.liftingYaxis,ScreenInfo.xmid + loc_pixel];
    num_dots = 3000;
    RNcoordinates = round(rand(num_dots,2) .* repmat([ScreenInfo.yaxis/replictation-1,ScreenInfo.xaxis/replictation-1],num_dots,1)) + 1; % the -1 and +1 was in order to prevent index = 0;
    
    coherPDF = mvnpdf(RNcoordinates,targetLoc./replictation,[40,0;0,80]);
    coherPDF = coherPDF ./ sum(coherPDF) .*80;
    num_frames = VSinfo.numFrames;
    frameBin = randi([0 1],num_dots,num_frames);
    wavelength = 60;
    oneRep = ones(wavelength,1);
    ideal = repmat(oneRep,num_frames./wavelength,1);
    include = 0;
    for f = 1:num_frames
        randInx = rand(num_dots,1)<coherPDF;
        frameBin(randInx,f) = ideal(f);
        frameBin(~randInx,f) = randi([0,1],sum(~randInx),1);
        include = include + sum(randInx);
    end
    
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