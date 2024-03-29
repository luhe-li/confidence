function Resp = LocalizeVisualStim(i, ExpInfo,...
    ScreenInfo,VSinfo,windowPtr)

%% precompute visual stimuli

%first compute the location of the visual stimulus in pixels
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
Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
vbl = Screen('Flip',windowPtr);
Screen('Flip',windowPtr, vbl + (VSinfo.numFrames - 0.5) * ScreenInfo.ifi);

% mask
for jj = 1:VSinfo.numFramesMasker 
Screen('DrawTexture', windowPtr, VSinfo.gwn_texture(rem(jj,VSinfo.GWNnumFrames)+1),[],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end
     
% % display visual stimulus
% for j = 1:VSinfo.numFrames 
%     Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
%         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%     Screen('Flip',windowPtr);
% end
% 
% % display post-stimulus mask
% for jj = 1:VSinfo.numFramesMasker
%     Screen('DrawTexture', windowPtr, VSinfo.gwn_texture(jj),[],...
%         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%     Screen('Flip',windowPtr);
% end
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);

% blank screen 2
% WaitSecs(ExpInfo.tBlank2);

%% response

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
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
% confidence response
Screen('TextSize',windowPtr,15);
SetMouse(x*2, yLoc*2, windowPtr);
HideCursor;
buttons = 0;
WaitSecs(0.2);
tic;
while sum(buttons)==0
    [conf_x,~,buttons] = GetMouse(windowPtr);
    conf_radius = abs(conf_x - x);
    potentialconfRcm = conf_radius/ScreenInfo.numPixels_perCM;
    potentialPoint = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * potentialconfRcm, ExpInfo.minPoint);
    
    potentialEnclosed = abs(ExpInfo.speakerLocCM(ExpInfo.randVisIdx(i)) - Resp.response_cm) <= potentialconfRcm;
    
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc+3, x, yLoc-3, 1);
    Screen('DrawLine', windowPtr, [255 255 255],x-conf_radius, yLoc, x+conf_radius, yLoc, 1);
    
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
Resp.conf_x = conf_x;
Resp.conf_radius_pixel= conf_radius;
Resp.conf_radius_cm  = Resp.conf_radius_pixel/ScreenInfo.numPixels_perCM;

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

% calculate points
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_cm = ExpInfo.speakerLocCM(Resp.target_idx);
Resp.target_pixel = Resp.target_cm * ScreenInfo.numPixels_perCM;
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