function Resp = LocalizeBothStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)

%% generate stimulus

% randomly draw VSinfo.num_randomDots [x;y] coordinates based on the centroid
% first compute the location of the visual stimulus in pixels
loc_pixel = round(ExpInfo.randVisPixel(i));
targetLoc = ScreenInfo.xmid + loc_pixel;
RNcoordinates = randn(2,1);
dots_targetLoc_coordinates = [targetLoc+(...
    ScreenInfo.numPixels_perCM.*VSinfo.SD_blob.*RNcoordinates(1,:));...
    ScreenInfo.liftingYaxis+(ScreenInfo.numPixels_perCM.*...
    VSinfo.SD_yaxis.*RNcoordinates(2,:))];
while 1
    %randomly draw 10 (x,y) coordinates based on the centroid
    RNcoordinates = randn(2,1);
    new_dot_targetLoc_coordinates = [targetLoc+(...
        ScreenInfo.numPixels_perCM.*VSinfo.SD_blob.*RNcoordinates(1,:));...
        ScreenInfo.liftingYaxis+(ScreenInfo.numPixels_perCM.*...
        VSinfo.SD_yaxis.*RNcoordinates(2,:))];
    dots_targetLoc_coordinates = [dots_targetLoc_coordinates,new_dot_targetLoc_coordinates];
    %make sure the center of the 10 blobs are aligned with the
    %predetermined location of the test stimulus
    dots_targetLoc_coordinates_shifted = shiftDotClouds(...
        dots_targetLoc_coordinates,loc_pixel,ScreenInfo);
    
    %check if they are within the boundaries
    check_withinTheLimit = CheckWithinTheBoundaries(...
        dots_targetLoc_coordinates_shifted,VSinfo.boxSize,ScreenInfo);
    
    %if the generated dots are within boundaries, then pass the
    %coordinates to the function generateDotClouds that gives out the
    %image texture.
    if check_withinTheLimit == 1
        if size(dots_targetLoc_coordinates,2) == VSinfo.num_randomDots
            dotClouds_targetLoc = generateDotClouds(windowPtr,...
                dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
            break;
        end
    else
        dots_targetLoc_coordinates = dots_targetLoc_coordinates(:,1:end-1);
    end
end
Resp.coordinates = dots_targetLoc_coordinates;
Resp.dot_x_mu_pixel =  mean(dots_targetLoc_coordinates(1,:));
Resp.dot_y_mu_pixel =  mean(dots_targetLoc_coordinates(2,:));
Resp.dot_x_sd_pixel = sqrt(sum((dots_targetLoc_coordinates(1,:) - Resp.dot_x_mu_pixel).^2)/size(dots_targetLoc_coordinates,2));
Resp.dot_y_sd_pixel = sqrt(sum((dots_targetLoc_coordinates(2,:) - Resp.dot_y_mu_pixel).^2)/size(dots_targetLoc_coordinates,2));

%% trial location info

% switch post-cue
if ExpInfo.randAVIdx(3,i) == 1 % A
    Resp.target_idx = ExpInfo.randAudIdx(i);
    Resp.target_cm = ExpInfo.randAudCM(i);
    Resp.target_deg = ExpInfo.randAudVA(i);
    Resp.target_pixel = ExpInfo.randAudPixel(i);
    Resp.post_cue = ExpInfo.randAVIdx(3,i);
else % V
    Resp.target_idx = ExpInfo.randVisIdx(i);
    Resp.target_cm = ExpInfo.randVisCM(i);
    Resp.target_deg = ExpInfo.randVisVA(i);
    Resp.target_pixel = ExpInfo.randVisPixel(i);
    Resp.post_cue = ExpInfo.randAVIdx(3,i);
end

%% trial start

% fixation
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tFixation);

% blank screen 1
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank1);

% present both stimulus
Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudIdx(i)),'>']; %arduino takes input in this format
fprintf(Arduino,input_on);
PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
vbl = Screen('Flip',windowPtr); % v onset
PsychPortAudio('Start',pahandle,1,0,0);
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5) * ExpInfo.tIFI);
WaitSecs(0.1);
input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudIdx(i)),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);
WaitSecs(0.1);

%% response

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
SetMouse(randi(ScreenInfo.xaxis,1), yLoc, windowPtr);
Screen('TextSize',windowPtr,15);
HideCursor;
buttons = 0;
tic;
while sum(buttons)==0
    [x,~,buttons] = GetMouse(windowPtr);
    HideCursor;
    x = min(x, ScreenInfo.xmid*2); 
    x = max(0,x);
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc-3, x, yLoc+3, 1);
    DrawFormattedText(windowPtr, ExpInfo.cue{ExpInfo.randAVIdx(3,i)}, 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-12,x+20,yLoc-6]);
    Screen('Flip',windowPtr);
    [~, ~, keyCode] = KbCheck(-1);
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    end
end
Resp.RT1  = toc;
Resp.response_pixel = x - ScreenInfo.xmid;
Resp.response_cm    = Resp.response_pixel/ScreenInfo.numPixels_perCM;
Resp.response_deg   = rad2deg(atan(Resp.response_cm/ExpInfo.sittingDistance));
HideCursor;

% confidence response
SetMouse(x*2, yLoc*2, windowPtr);
HideCursor;
WaitSecs(0.2);
buttons = 0;
tic;
pm = PsychPowerMate('Open');
[~, dialPos] = PsychPowerMate('Get',pm); %initalize powermate
initDialPos = dialPos;
while sum(buttons)==0

    [~,~,buttons] = GetMouse(windowPtr);
    [~, dialPos] = PsychPowerMate('Get',pm); %update dial postion
    
    conf_radius = ExpInfo.dialScaler * abs(dialPos - initDialPos);
    potentialconfRcm = conf_radius/ScreenInfo.numPixels_perCM;
    potentialPoint = max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * potentialconfRcm, ExpInfo.minPoint);
    potentialEnclosed = abs(Resp.target_cm - Resp.response_cm) <= potentialconfRcm;
   
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
    Screen('FillRect', windowPtr, [255 255 255]./7, [x-conf_radius, yLoc-ExpInfo.conf_bar_height/2 + 3, x+conf_radius, yLoc+ExpInfo.conf_bar_height/2 - 3]);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc+ExpInfo.conf_bar_height/2, x-conf_radius, yLoc-ExpInfo.conf_bar_height/2, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x+conf_radius, yLoc+ExpInfo.conf_bar_height/2, x+conf_radius, yLoc-ExpInfo.conf_bar_height/2, 1);   
    
    DrawFormattedText(windowPtr, ['Potential score: ' num2str(round(potentialPoint,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-12,x+20,yLoc-6]);
    
    Screen('Flip',windowPtr);
    [~, ~, keyCode] = KbCheck(-1);
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    end
end
Resp.RT2             = toc;
Resp.conf_radius_pixel= conf_radius;
Resp.conf_radius_cm  = Resp.conf_radius_pixel/ScreenInfo.numPixels_perCM;

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
if ~rem(i,3) && ~ismember(i-1, ExpInfo.breakTrials) && ~ismember(i-2, ExpInfo.breakTrials) || ExpInfo.practice % every 3 trials give feedback OR it's practice
    Screen('TextSize',windowPtr,25);
    DrawFormattedText(windowPtr, ['Score of the last trial: ' num2str(round(potentialPoint * potentialEnclosed,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [ScreenInfo.xmid-20,yLoc-3,ScreenInfo.xmid+20,yLoc+3]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.tFeedback);
end
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tITI);

%calculate points
Resp.enclosed = abs(Resp.target_pixel - Resp.response_pixel) <= Resp.conf_radius_pixel;
bestRadius_cm = abs(Resp.target_cm - Resp.response_cm);
Resp.maxPtPossible = max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * bestRadius_cm, ExpInfo.minPoint);
if Resp.enclosed
    Resp.point = Resp.maxPtPossible;
else
    Resp.point = 0;
end

end