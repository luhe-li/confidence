function Resp = LocalizeVisualStim(i, ExpInfo, ScreenInfo,VSinfo,windowPtr)

%% make visual stimuli

blob_coordinates = [ExpInfo.randVisPixel(i)+ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

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
for jj = 1:ExpInfo.tStimFrame
Screen('DrawTexture', windowPtr, dotCloud, [],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);

%% collect response

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
SetMouse(randi(ScreenInfo.xaxis,1), yLoc, windowPtr);
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

% confidence response
Screen('TextSize',windowPtr,15);
SetMouse(x*2, yLoc*2, windowPtr);
HideCursor;
WaitSecs(0.2);
tic;
pm = PsychPowerMate('Open');
[buttonPM, dialPos] = PsychPowerMate('Get',pm); %initalize powermate
initDialPos = dialPos;
while ~buttonPM
    [buttonPM, dialPos] = PsychPowerMate('Get',pm); %update dial postion
    
    conf_radius = dialScaler * abs(dialPos - initDialPos);
    potentialconfRcm = conf_radius/ScreenInfo.numPixels_perCM;
    potentialPoint = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * potentialconfRcm, ExpInfo.minPoint);
    
    potentialEnclosed = abs(ExpInfo.speakerLocCM(ExpInfo.randVisIdx(i)) - Resp.response_cm) <= potentialconfRcm;
    
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
    Screen('FillRect', windowPtr, [255 255 255]./7, [x-conf_radius, yLoc-height/2 + 3, x+conf_radius, yLoc+height/2 - 3]);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc+height/2, x-conf_radius, yLoc-height/2, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x+conf_radius, yLoc+height/2, x+conf_radius, yLoc-height/2, 1);   
    
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
if ~rem(i,3) || ExpInfo.practice == 2 % every three trials give feedback, if not practice
    DrawFormattedText(windowPtr, ['Score of the last trial: ' num2str(round(potentialPoint * potentialEnclosed,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [ScreenInfo.xmid-20,yLoc-3,ScreenInfo.xmid+20,yLoc+3]);
    WaitSecs(0.1);
end
Screen('Flip',windowPtr);

% calculate points
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_pixel = ExpInfo.randVisPixel(i);
Resp.target_cm = ExpInfo.randVisCM(i);
Resp.target_deg = ExpInfo.randVisVA(i);
Resp.enclosed = abs(Resp.target_pixel - Resp.response_pixel) <= Resp.conf_radius_cm;
bestRadius_pixel = abs(Resp.target_pixel - Resp.response_pixel);
Resp.maxPtPossible = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * bestRadius_pixel, ExpInfo.minPoint);
if Resp.enclosed
    Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.bestRadius_pixel, ExpInfo.minPoint);
else
    Resp.point = 0;
end

end