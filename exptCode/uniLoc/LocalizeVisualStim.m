function Resp = LocalizeVisualStim(i, ExpInfo,...
    ScreenInfo,VSinfo,windowPtr)

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

%%  present visual stimulus

%first compute the location of the visual stimulus in pixels
loc_pixel = round(ExpInfo.loc_pixel(i));
targetLoc = ScreenInfo.xmid + loc_pixel;

while 1
    %randomly draw 10 (x,y) coordinates based on the centroid
    RNcoordinates = randn(2,VSinfo.num_randomDots);
    dots_targetLoc_coordinates = vertcat(targetLoc+(...
        ScreenInfo.numPixels_perCM.*VSinfo.SD_blob.*RNcoordinates(1,:)),...
        ScreenInfo.liftingYaxis+(ScreenInfo.numPixels_perCM.*...
        VSinfo.SD_yaxis.*RNcoordinates(2,:)));

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
        dotClouds_targetLoc = generateDotClouds(windowPtr,...
            dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
        break;
    end
end


% display visual stimulus

for j = 1:VSinfo.numFrames %100 ms
    Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
end

% blank screen 2
WaitSecs(ExpInfo.tBlank2);

%% response

% perceptual response
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


% confidence response
SetMouse(x, ScreenInfo.ymid*2, windowPtr);
buttons = 0;
WaitSecs(0.2);
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
Resp.loc_idx = ExpInfo.loc_idx(i);
Resp.loc_cm  = ExpInfo.loc_cm(i);
Resp.loc_deg = ExpInfo.loc_deg(i);
Resp.loc_pixel = ExpInfo.loc_pixel(i);
Resp.enclosed = Resp.loc_cm >= Resp.response_cm - Resp.conf_radius_cm & ...
    Resp.loc_cm <= Resp.response_cm + Resp.conf_radius_cm;
if Resp.enclosed
    Resp.point = 0.01 * max(ExpInfo.maxPoint - ExpInfo.dropRate * 2 * Resp.conf_radius_cm, ExpInfo.minPoint);
else
    Resp.point = 0;
end

end