function Resp = LocalizeAuditoryStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)
   
    % fixation
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
        ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
        ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.tFixation);

    % blank screen 1
    Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.tBlank1);

    % present stimulus
    input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudLoc(i)),'>']; %arduino takes input in this format
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(ExpInfo.tStim);
    input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudLoc(i)),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);

    % blank screen 2
%     Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
%         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%     Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.tBlank2);

    % perception response
    yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
    SetMouse(randi(ScreenInfo.xmid*2,1), ScreenInfo.ymid*2, windowPtr);
    buttons = 0;
    tic
    while sum(buttons)==0
        [x,~,buttons] = GetMouse(windowPtr); HideCursor;
        x = min(x, ScreenInfo.xmid*2); x = max(0,x);
        Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
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
    Response_pixel = x;
    Resp.Response_cm    = (Response_pixel -  ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
    Resp.Response_deg   = rad2deg(atan(Resp.Response_cm/ExpInfo.sittingDistance));
    
    
    % wo xie de bu fen
    SetMouse(x, ScreenInfo.ymid*2, windowPtr);
    buttons = 0;
    WaitSecs(0.2)
    tic
    while sum(buttons)==0
        [conf_x,~,buttons] = GetMouse(windowPtr); HideCursor;
        Resp.conf_radius = abs(conf_x - x);
        
        Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
            [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
        Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
        Screen('DrawLine', windowPtr, [255 255 255],x-Resp.conf_radius, yLoc, x+Resp.conf_radius, yLoc, 1);
        Screen('Flip',windowPtr);
        [~, ~, keyCode] = KbCheck();
        if keyCode(KbName('ESCAPE'))
            sca;
            error('Escape');
        end
    end
    Resp.RT2             = toc;
    
    % confidence response

    % ITI
    Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    WaitSecs(ExpInfo.ITI);

end