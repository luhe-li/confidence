function Resp = RndBinVisStim(i, ExpInfo, ScreenInfo,VSinfo,windowPtr)

%% precompute visual stimuli
height = 200;
noise_sd = 20;
stim_sd = VSinfo.SD_blob(i) .* 6;
wBack = 0.5;
stimTx = generateRippleStim(VSinfo,ExpInfo,ScreenInfo,windowPtr,i, height, noise_sd, stim_sd, wBack);


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

% display visual stimulus
for jj = 1:VSinfo.numFrames
Screen('DrawTexture', windowPtr, stimTx(jj),[],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);

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

% Record target location
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_cm = ExpInfo.randVisCM(i);
Resp.target_pixel = ExpInfo.randVisPixel(i);
Resp.target_deg = ExpInfo.randVisVA(i);

end