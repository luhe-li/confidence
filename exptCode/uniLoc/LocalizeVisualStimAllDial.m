function Resp = LocalizeVisualStimAllDial(i, ExpInfo, ScreenInfo,VSinfo,windowPtr)

%% precompute visual stimuli
height = 200;
noise_sd = 20;
stim_sd = VSinfo.SD_blob(i) .* 8;
wBack = 0.55;
stimTx = generateRippleStim(VSinfo,ExpInfo,ScreenInfo,windowPtr,i, height, noise_sd, stim_sd, wBack);

dialScalerConf = 2;
dialScalerEst = 10;
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
pm = PsychPowerMate('Open');
[buttonPM, ~] = PsychPowerMate('Get',pm); %initalize powermate
initDialPos = -randi(ceil(ScreenInfo.xaxis/dialScalerEst),1);

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
HideCursor;
tic;
while ~buttonPM
    [buttonPM, dialPos] = PsychPowerMate('Get',pm); %update dial postion
    x = dialScalerEst * abs(dialPos - initDialPos);
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

[buttonPM, dialPos] = PsychPowerMate('Get',pm); %initalize powermate
initDialPos = dialPos;
while ~buttonPM
    [buttonPM, dialPos] = PsychPowerMate('Get',pm); %update dial postion
    
    conf_radius = dialScalerConf * abs(dialPos - initDialPos);
    potentialconfRcm = conf_radius/ScreenInfo.numPixels_perCM;
    potentialconfRdeg = rad2deg(atan(potentialconfRcm/ExpInfo.sittingDistance));
    potentialPoint = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * potentialconfRdeg, ExpInfo.minPoint);
    
    potentialEnclosed = abs(ExpInfo.speakerLocVA(ExpInfo.randVisIdx(i)) - Resp.response_deg) <= potentialconfRdeg;
    
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc+height/2, x-conf_radius, yLoc-height/2, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x+conf_radius, yLoc+height/2, x+conf_radius, yLoc-height/2, 1);
    
    if ExpInfo.practice == 2
    DrawFormattedText(windowPtr, ['Actual score: ' num2str(round(potentialPoint * potentialEnclosed,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-25,x+20,yLoc-19]);
    DrawFormattedText(windowPtr, ['Potential score: ' num2str(round(potentialPoint,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-12,x+20,yLoc-6]);
    else
    DrawFormattedText(windowPtr, num2str(round(potentialPoint,2)), 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [x-20,yLoc-12,x+20,yLoc-6]);
    end
    
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
if ~rem(i,3) && ExpInfo.practice ~= 2 % every three trials give feedback, if not practice
    DrawFormattedText(windowPtr, ['Score: ' num2str(round(potentialPoint * potentialEnclosed,2))], 'center', 'center', ...
        [255 255 255],[], [], [], [], [], ...
        [ScreenInfo.xmid-20,yLoc-3,ScreenInfo.xmid+20,yLoc+3]);
end
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

% calculate points
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_cm = ExpInfo.speakerLocCM(Resp.target_idx);
Resp.target_deg = rad2deg(atan(Resp.target_cm/ExpInfo.sittingDistance));
Resp.enclosed = abs(Resp.target_cm - Resp.response_cm) <= Resp.conf_radius_cm;
bestRadius_cm = abs(Resp.target_cm - Resp.response_cm);
Resp.maxPtPossible = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * bestRadius_cm, ExpInfo.minPoint);
if Resp.enclosed
    Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
else
    Resp.point = 0;
end

end