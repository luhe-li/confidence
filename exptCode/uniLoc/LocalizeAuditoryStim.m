function Resp = LocalizeAuditoryStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)
   
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
    input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudLoc(i)),'>']; %arduino takes input in this format
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(ExpInfo.tStim);
    input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudLoc(i)),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);

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
    
    disp(x)
    % confidence response
    SetMouse(x, ScreenInfo.ymid, windowPtr);
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
    Resp.enclosed = Resp.loc_cm >= Resp.response_cm - Resp.conf_radius_cm & ...
        Resp.loc_cm <= Resp.response_cm + Resp.conf_radius_cm;
    if Resp.enclosed
        Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
    else
        Resp.point = 0;
    end
        
end