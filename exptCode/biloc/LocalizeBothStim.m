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
Screen('Flip',windowPtr, vbl + (ExpInfo.frameStim - 0.5) * ScreenInfo.ifi);
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

% blank screen 2
% WaitSecs(ExpInfo.tBlank2);

%% response

% perception response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
Screen('TextSize',windowPtr,20);
SetMouse(randi(ScreenInfo.xmid*2,1), yLoc*2, windowPtr);
buttons = 0;
tic;
HideCursor;
while sum(buttons)==0
    [x,~,buttons] = GetMouse(windowPtr); 
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
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'Are you confident about your estimation?\nYes: 1\nNo: 2', ...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,[255 255 255]);
Screen('Flip',windowPtr);

resp=1; tic;
while resp
    [~, ~, keyCode] = KbCheck();
    if keyCode(KbName('numLock'))
        ShowCursor;
        sca;
        error('Escape');
    elseif keyCode(KbName('1'))
        conf = 1;
        resp = 0;
    elseif keyCode(KbName('2'))
        conf = 0;
        resp = 0;
    end
end
Resp.RT2 = toc;
Resp.conf = conf;

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
% Resp.enclosed = abs(Resp.target_cm - Resp.response_cm) <= Resp.conf_radius_cm;
% bestRadius_cm = abs(Resp.target_cm - Resp.response_cm);
% Resp.maxPtPossible = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * bestRadius_cm, ExpInfo.minPoint);
% if Resp.enclosed
%     Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
% else
%     Resp.point = 0;
% end

end