function Resp = LocalizePointStim(i, ExpInfo,...
    ScreenInfo,VSinfo,windowPtr)

%% precompute visual stimuli

%first compute the location of the visual stimulus in pixels
jitter = VSinfo.jitter_lb + (VSinfo.jitter_ub - VSinfo.jitter_lb) * rand;
loc_pixel = round(ExpInfo.randVisPixel(i) + jitter);
xLoc = ScreenInfo.xmid + loc_pixel;
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
targetLoc = [xLoc, yLoc];
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

% display visual stimulu
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('DrawDots',windowPtr,targetLoc,5,[255 255 255],[],1);
vbl = Screen('Flip',windowPtr);
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr, vbl + (VSinfo.numFrames - 0.5) * ScreenInfo.ifi);

% blank screen 2
WaitSecs(ExpInfo.tBlank2);

%% response

% perceptual response

SetMouse(randi(ScreenInfo.xaxis*2,1), yLoc*2, windowPtr);
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
HideCursor;

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

% calculate points
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_pixel = loc_pixel;
Resp.target_cm    = loc_pixel./ScreenInfo.numPixels_perCM;
Resp.target_deg = rad2deg(atan(Resp.target_cm/ExpInfo.sittingDistance));

end