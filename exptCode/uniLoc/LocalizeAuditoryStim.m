function Resp = LocalizeAuditoryStim(i, ExpInfo,...
    ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr)

%% start the trial

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
input_on = ['<',num2str(1),':',num2str(ExpInfo.randAudIdx(i)),'>']; %arduino takes input in this format
fprintf(Arduino,input_on);
PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
PsychPortAudio('Start',pahandle,1,0,0);
switch ExpInfo.mode
    case 1 % experiment mode
        WaitSecs(0.1);
    case 2 % debug mode
        WaitSecs(3.5);
end
input_off = ['<',num2str(0),':',num2str(ExpInfo.randAudIdx(i)),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);
WaitSecs(0.1);

%% response

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
SetMouse(randi(ScreenInfo.xaxis*2,1), yLoc*2, windowPtr);
HideCursor;
resp = 1;
tic;
stopRecorded = 0;
x = -1;
while resp
    cache = x;
    [x,~,~] = GetMouse(windowPtr);
    HideCursor;
    x = min(x, ScreenInfo.xmid*2);
    x = max(0,x);
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc-3, x, yLoc+3, 1);
    Screen('Flip', windowPtr);

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

Resp.target_idx = ExpInfo.randAudIdx(i); % visual location that corresponds to speaker index
Resp.target_cm = ExpInfo.randAudCM(i);
Resp.target_pixel = ExpInfo.randAudPixel(i);
Resp.target_deg = ExpInfo.randAudVA(i);

end