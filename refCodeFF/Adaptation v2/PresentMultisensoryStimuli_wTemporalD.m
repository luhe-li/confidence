function [localization, RT1, unity, RT2] = PresentMultisensoryStimuli_wTemporalD(...
    TrialNum,ExpInfo,ScreenInfo,VSinfo, AudInfo,motorArduino,numRailSteps,...
    pahandle,windowPtr)
    %----------------------------------------------------------------------
    %-----------Calculate the coordinates of the target stimuli------------
    %----------------------------------------------------------------------
    %display visual stimuli
    targetLoc = ScreenInfo.xmid + ScreenInfo.numPixels_perCM.*...
                VSinfo.arrangedLocs_cm(TrialNum);
    
    %Make visual stimuli
    blob_coordinates = [targetLoc, ScreenInfo.liftingYaxis];    
    dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

    %----------------------------------------------------------------------
    %--------------Move the motor to the correct location------------------
    %----------------------------------------------------------------------
    %display Mask Noise
    PsychPortAudio('FillBuffer', pahandle, AudInfo.MaskNoise);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    
    %calculate the wait time
    movingSteps = AudInfo.moving_locations_steps(TrialNum,:); 
    waitTime1   = FindWaitTime(movingSteps);      
    %move the speaker to the location we want 
    for s = 1:length(movingSteps)
        if movingSteps(s) < 0  %when AuditoryLoc is negative, move to the left
            fprintf(motorArduino,['%c','%d'], ['p', numRailSteps*abs(movingSteps(s))]);
        else
            fprintf(motorArduino,['%c','%d'], ['n', numRailSteps*movingSteps(s)]);
        end
        %wait shortly
        WaitSecs(waitTime1(s));
    end
    %wait shortly
    WaitSecs(AudInfo.waitTime);
    PsychPortAudio('Stop', pahandle);

    %----------------------------------------------------------------------
    %---------------------display audiovisual stimuli----------------------
    %----------------------------------------------------------------------
    %show fixation cross for 1 s and then a blank screen for 2 s
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
        ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
        ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
    Screen('Flip',windowPtr); WaitSecs(0.1);
    Screen('Flip',windowPtr); WaitSecs(0.5);

    if VSinfo.timing_relative(TrialNum) < 0 %present V first 
        for j = 1:VSinfo.numFrames 
            Screen('DrawTexture',windowPtr,dotCloud,[],[0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            Screen('Flip',windowPtr);
        end 
        %temporal offset = 250ms
        Screen('Flip',windowPtr); WaitSecs(abs(VSinfo.timing_relative(TrialNum))/1e3);
        %start playing the sound while displaying the visual stimulus
        PsychPortAudio('FillBuffer', pahandle, AudInfo.GaussianWhiteNoise);
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        WaitSecs(AudInfo.adaptationDuration);
        PsychPortAudio('Stop', pahandle);
        %black screen for 250ms
        Screen('Flip',windowPtr); WaitSecs(abs(VSinfo.timing_relative(TrialNum))/1e3);
    else
        %start playing the sound while displaying the visual stimulus
        PsychPortAudio('FillBuffer', pahandle, AudInfo.GaussianWhiteNoise);
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        WaitSecs(AudInfo.adaptationDuration);
        PsychPortAudio('Stop', pahandle);
        %black screen for 250ms
        Screen('Flip',windowPtr); WaitSecs(abs(VSinfo.timing_relative(TrialNum))/1e3);
        for j = 1:VSinfo.numFrames 
            Screen('DrawTexture',windowPtr,dotCloud,[],[0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            Screen('Flip',windowPtr);
        end 
    end
    
    %----------------------------------------------------------------------
    %--------------Record response using the pointing device---------------
    %----------------------------------------------------------------------
    Screen('Flip',windowPtr); WaitSecs(0.5);
    yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
    Screen('TextSize', windowPtr, 12);
    modality = ExpInfo.localize_modality(TrialNum);
    SetMouse(randi(ScreenInfo.xmid*2,1), yLoc, windowPtr); 
    buttons = zeros(1,16); tic
    %localize the stimulus using a visual cursor
    while buttons(1) == 0
        [x,~,buttons] = GetMouse; HideCursor;
        Screen('FillRect', windowPtr, ScreenInfo.cursorColor(modality,:),...
            [x-3 yLoc-24 x+3 yLoc-12]);
        Screen('FillPoly', windowPtr, ScreenInfo.cursorColor(modality,:),...
            [x-3 yLoc-12; x yLoc; x+3 yLoc-12]);
        Screen('DrawText', windowPtr, ScreenInfo.dispModality(modality),...
            x-5, yLoc-30,[255 255 255]);
        Screen('Flip',windowPtr);
    end
    RT1            = toc;
    Response_pixel = x;
    Response_cm    = (Response_pixel - ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
    localization   = rad2deg(atan(Response_cm/ExpInfo.sittingDistance));
    Screen('Flip',windowPtr); WaitSecs(0.1);
    
    %Unity judgment
    if ExpInfo.bool_unityReport(TrialNum) == 1
        Screen('TextSize', windowPtr, 25);
        click = 0; tic
        while click == 0
            Screen('FillRect', windowPtr, [0,0,0],[x-3 yLoc-100 x+3 yLoc-100]);
            Screen('DrawText', windowPtr, 'C=1    OR    C=2',ScreenInfo.xmid-87.5,...
                yLoc-10,[255 255 255]);
            Screen('Flip',windowPtr);
            %click left button: C=1; click right button: C=2
            [click,~,~,unity] = GetClicks;
        end
        RT2     = toc;
        %show a white frame to confirm the choice
        x_shift = ScreenInfo.x_box_unity(unity,:);
        Screen('DrawText', windowPtr, 'C=1    OR    C=2',ScreenInfo.xmid-87.5,...
                yLoc-10,[255 255 255]);
        Screen('FrameRect', windowPtr, [255,255,255], [ScreenInfo.xmid+x_shift(1),...
            yLoc-10+ScreenInfo.y_box_unity(1), ScreenInfo.xmid+x_shift(2),...
            yLoc-10+ScreenInfo.y_box_unity(2)]);
        Screen('Flip',windowPtr); WaitSecs(0.1);
        Screen('Flip',windowPtr); WaitSecs(0.1);
    else
        unity = NaN; RT2 = NaN; 
    end
end
    
    