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
input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudIdx(i)),'>']; %arduino takes input in this format
fprintf(Arduino,input_on);
PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
vbl = Screen('Flip',windowPtr); % v onset
PsychPortAudio('Start',pahandle,0,0,0);
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr, vbl + (ExpInfo.frameStim - 0.5) * ScreenInfo.ifi);
input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudIdx(i)),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);

% blank screen 2
WaitSecs(ExpInfo.tBlank2);

%% response

% perception response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
Screen('TextSize',windowPtr,14);
SetMouse(randi(ScreenInfo.xmid*4,1), yLoc*2, windowPtr);
buttons = 0;
tic;
HideCursor;
while sum(buttons)==0
    [x,~,buttons] = GetMouse(windowPtr); 
    x = min(x, ScreenInfo.xmid*2); x = max(0,x);
    HideCursor;
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
Resp.response_pixel = x;
Resp.response_cm    = (Resp.response_pixel -  ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
Resp.response_deg   = rad2deg(atan(Resp.response_cm/ExpInfo.sittingDistance));

% confidence response
SetMouse(x*2, yLoc*2, windowPtr);
buttons = 0;
WaitSecs(0.2);
tic;
while sum(buttons)==0
    [conf_x,~,buttons] = GetMouse(windowPtr); 
    conf_radius = abs(conf_x - x);
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
    DrawFormattedText(windowPtr,ExpInfo.cue{ExpInfo.randAVIdx(3,i)}, 'center', 'center', ...
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
Resp.RT2  = toc;
Resp.conf_radius_pixel= conf_radius;
Resp.conf_radius_cm  = Resp.conf_radius_pixel/ScreenInfo.numPixels_perCM;

% Common Cause response
Screen('TextSize',windowPtr,30);
ccSliderLength = 200; % common cause slide length
ccSliderHeight = 26; % common cause slide length
ccSliderRect([1,3]) = ScreenInfo.xmid-1 + [-ccSliderLength/2+1,ccSliderLength/2];
ccSliderRect([2,4]) = yLoc + [-ccSliderHeight/2,ccSliderHeight/2];
SetMouse(randi(ccSliderRect([1,3]),1)*2, yLoc*2, windowPtr);

buttons = 0;
WaitSecs(0.2);
tic;
while sum(buttons)==0
    [caus_x,~,buttons] = GetMouse(windowPtr);
    caus_x = min(caus_x, ccSliderRect(3));
    caus_x = max(ccSliderRect(1),caus_x);
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],caus_x, ccSliderRect(2), caus_x, ccSliderRect(4), 1);
    Screen('FrameRect', windowPtr, [255 255 255], ccSliderRect', 1);
    DrawFormattedText(windowPtr, 'Common', 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        ccSliderRect - [ccSliderLength*1.2,0,ccSliderLength*1.2,0]);
    DrawFormattedText(windowPtr, 'Separate', 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        ccSliderRect + [ccSliderLength*1.2,0,ccSliderLength*1.2,0]);
    Screen('Flip',windowPtr);
    [~, ~, keyCode] = KbCheck(-1);
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    end
end
Resp.RT3  = toc;
Resp.unityConf = (caus_x - ccSliderRect(1) ) / ccSliderLength * 2 -1;
Resp.unity = Resp.unityConf > 0;

%ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

%calculate points
if ExpInfo.randAVIdx(3,i) == 1 % A
    Resp.target_idx = ExpInfo.randAVIdx(1,i);
else % V
    Resp.target_idx = ExpInfo.randAVIdx(2,i);
end
Resp.target_cm = ExpInfo.speakerLocCM(Resp.target_idx);
Resp.target_deg = rad2deg(atan(Resp.target_cm/ExpInfo.sittingDistance));
Resp.enclosed = abs(Resp.target_cm - Resp.response_cm) <= Resp.conf_radius_cm;
if Resp.enclosed
    Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
else
    Resp.point = 0;
end

end