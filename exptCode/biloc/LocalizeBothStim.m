function Resp = LocalizeBothStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)

%% randomly draw VSinfo.num_randomDots [x;y] coordinates based on the centroid
% first compute the location of the visual stimulus in pixels
loc_pixel = round(ExpInfo.loc_pixel(i));
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
        else
            dots_targetLoc_coordinates = [dots_targetLoc_coordinates,new_dot_targetLoc_coordinates];
        end
    end
end
%%

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

    % present auditory stimulus
    Resp.visFrameTime = NaN(1,VSinfo.numFrames);
    input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudLoc(i)),'>']; %arduino takes input in this format
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    Resp.audStartTime = PsychPortAudio('Start',pahandle,0,0,1,GetSecs+AudInfo.stimDura);
    % WaitSecs(ExpInfo.tStim);
    for j = 1:VSinfo.numFrames %100 ms, 6 frames
        Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Resp.visFrameTime(j) = Screen('Flip',windowPtr);
    end
    input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudLoc(i)),'>'];
    fprintf(Arduino,input_off);
    % PsychPortAudio('Stop',pahandle);

    % blank screen 2
    WaitSecs(ExpInfo.tBlank2);

    % perception response
    yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
    SetMouse(randi(ScreenInfo.xmid*2,1), ScreenInfo.ymid*2, windowPtr);
    buttons = 0;
    tic
    while sum(buttons)==0
        [x,~,buttons] = GetMouse(windowPtr); HideCursor;
        x = min(x, ScreenInfo.xmid*2); x = max(0,x);
        Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
                [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
        Screen('DrawLine', windowPtr, [255 255 255],x, yLoc-3, x, yLoc+3, 1);
        Screen('Flip',windowPtr);
        [~, ~, keyCode] = KbCheck();
        if keyCode(KbName('ESCAPE'))
           sca;
           error('Escape');
        end
    end
    Resp.RT1  = toc;
    HideCursor;
    Resp.response_pixel = x;
    Resp.response_cm    = (Resp.response_pixel -  ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
    Resp.response_deg   = rad2deg(atan(Resp.response_cm/ExpInfo.sittingDistance));
    
    
    % confidence response
    SetMouse(x, ScreenInfo.ymid*2, windowPtr);
    buttons = 0;
    WaitSecs(0.2)
    tic
    while sum(buttons)==0
        [conf_x,~,buttons] = GetMouse(windowPtr); HideCursor;
        conf_radius = abs(conf_x - x);
        
        Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
            [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
        Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
        Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
        Screen('Flip',windowPtr);
        [~, ~, keyCode] = KbCheck();
        if keyCode(KbName('ESCAPE'))
            sca;
            error('Escape');
        end
    end
    Resp.RT2             = toc;
    HideCursor;
    Resp.conf_radius_pixel= conf_radius;
    Resp.conf_radius_cm  = Resp.conf_radius_pixel/ScreenInfo.numPixels_perCM;
  
    % ITI
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ITI);
  
    % calculate points
    Resp.loc_idx = ExpInfo.randAudLoc(i);
    Resp.loc_cm  = ExpInfo.loc_cm(i);
    Resp.loc_deg = rad2deg(atan(Resp.loc_cm/ExpInfo.sittingDistance));
    Resp.enclosed = abs(Resp.loc_cm - Resp.response_cm) <= Resp.conf_radius_cm;
    if Resp.enclosed
        Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
    else
        Resp.point = 0;
    end
        
end