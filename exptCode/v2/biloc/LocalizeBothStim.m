function Resp = LocalizeBothStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)

%% generate stimulus

% randomly draw VSinfo.num_randomDots [x;y] coordinates based on the centroid
% first compute the location of the visual stimulus in pixels
loc_pixel = round(ExpInfo.randVisPixel(i));
targetLoc = ScreenInfo.xmid + loc_pixel;
RNcoordinates = randn(2,1);
dots_targetLoc_coordinates = [targetLoc+(...
    ScreenInfo.numPixels_perCM.*VSinfo.SD_blob(i).*RNcoordinates(1,:));...
    ScreenInfo.liftingYaxis+(ScreenInfo.numPixels_perCM.*...
    VSinfo.SD_yaxis.*RNcoordinates(2,:))];
while 1
    %randomly draw 10 (x,y) coordinates based on the centroid
    RNcoordinates = randn(2,1);
    new_dot_targetLoc_coordinates = [targetLoc+(...
        ScreenInfo.numPixels_perCM.*VSinfo.SD_blob(i).*RNcoordinates(1,:));...
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
Resp.vStimDotsCoor = dots_targetLoc_coordinates;

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
PsychPortAudio('Start',pahandle,1,0,0);
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5) * ExpInfo.tIFI);
WaitSecs(0.1);
input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudIdx(i)),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);

% mask
for jj = 1:VSinfo.numFramesMasker 
Screen('DrawTexture', windowPtr, VSinfo.gwn_texture(rem(jj,VSinfo.GWNnumFrames)+1),[],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

%% response

% perception response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
Screen('TextSize',windowPtr,20);
SetMouse(randi(ScreenInfo.xaxis*2,1), yLoc*2, windowPtr);
HideCursor;
resp = 1;
tic;
stopRecorded = 0;
x = -1;
while resp
    cache = x;
    [x,~,~] = GetMouse(windowPtr); 
    x = min(x, ScreenInfo.xmid*2); 
    x = max(0,x);
    HideCursor;
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc-3, x, yLoc+3, 1);
    DrawFormattedText(windowPtr, ExpInfo.cue{ExpInfo.randAVIdx(3,i)}, 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-12,x+20,yLoc-6]);
    Screen('Flip',windowPtr);
    
    locdiff = abs(cache - x);
    if locdiff ~= 0
        stopRecorded = 0;
    elseif locdiff == 0 && ~stopRecorded
        mouseStopT = GetSecs();
        stopRecorded = 1;
    end
    
    % Check the keyboard
    [keyIsDown, startTime, keyCode] = KbCheck();
    if keyIsDown
        % Check if any of the specified keys are pressed
        if keyCode(KbName('a')) || keyCode(KbName('s')) || keyCode(KbName('d')) || keyCode(KbName('f'))
            [releaseTime, ~, ~] = KbReleaseWait();
            Resp.PressDuration = releaseTime - startTime;
            Resp.mouseStopDuration = startTime - mouseStopT;
            if keyCode(KbName('a'))
                conf = 1;
            elseif keyCode(KbName('s'))
                conf = 2;
            elseif keyCode(KbName('d'))
                conf = 3;
            elseif keyCode(KbName('f'))
                conf = 4;
            end
            resp = 0;
        end
        if keyCode(KbName('ESCAPE'))
            sca;
            ShowCursor;
            Screen('CloseAll');
            error('Escape!');
        end
    end
end
Resp.conf = conf;
Resp.RT1  = toc;
Resp.response_pixel = x;
Resp.response_cm    = (Resp.response_pixel -  ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
Resp.response_deg   = rad2deg(atan(Resp.response_cm/ExpInfo.sittingDistance));
HideCursor;

% ITI
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
Resp.target_pixel = Resp.target_cm .* ScreenInfo.numPixels_perCM;
Resp.target_deg = rad2deg(atan(Resp.target_cm/ExpInfo.sittingDistance));

end