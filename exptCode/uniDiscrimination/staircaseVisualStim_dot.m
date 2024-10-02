function Resp = staircaseVisualStim_dot(i, j, k, ExpInfo, ScreenInfo, VSinfo, windowPtr, Resp, flag_easy)

if ~exist('flag_easy','Var'); flag_easy = 0; end
if i == ExpInfo.n_trial; flag_easy = 1; end

%% precompute

% standard stimulus
x_loc_pixel = round(ExpInfo.standard_loc * ExpInfo.px_per_cm);
targetLoc = ScreenInfo.xmid + x_loc_pixel;
RNcoordinates = randn(2,1);
dots_targetLoc_coordinates = [targetLoc+(...
    ScreenInfo.px_per_cm.*VSinfo.SD_blob.*RNcoordinates(1,:));...
    ScreenInfo.liftingYaxis+(ScreenInfo.px_per_cm.*...
    VSinfo.SD_yaxis.*RNcoordinates(2,:))];
while 1
    %randomly draw 10 (x,y) coordinates based on the centroid
    RNcoordinates = randn(2,1);
    new_dot_targetLoc_coordinates = [targetLoc+(...
        ScreenInfo.px_per_cm.*VSinfo.SD_blob.*RNcoordinates(1,:));...
        ScreenInfo.liftingYaxis+(ScreenInfo.px_per_cm.*...
        VSinfo.SD_yaxis.*RNcoordinates(2,:))];
    dots_targetLoc_coordinates = [dots_targetLoc_coordinates,new_dot_targetLoc_coordinates];
    %make sure the center of the 10 blobs are aligned with the
    %predetermined location of the test stimulus
    dots_targetLoc_coordinates_shifted = shiftDotClouds(...
        dots_targetLoc_coordinates,x_loc_pixel,ScreenInfo);
    
    %check if they are within the boundaries
    check_withinTheLimit = CheckWithinTheBoundaries(...
        dots_targetLoc_coordinates_shifted,VSinfo.boxSize,ScreenInfo);
    
    %if the generated dots are within boundaries, then pass the
    %coordinates to the function generateDotClouds that gives out the
    %image texture.
    if check_withinTheLimit == 1
        if size(dots_targetLoc_coordinates,2) == VSinfo.num_randomDots
            std_stim = generateDotClouds(windowPtr,...
                dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
            break;
        end
    else
        dots_targetLoc_coordinates = dots_targetLoc_coordinates(:,1:end-1);
    end
end
Resp.standard_loc(k,j,i) = ExpInfo.standard_loc;

% comparison stimulus
x_loc_pixel = round(Resp.comparison_loc(k,j,i) * ExpInfo.px_per_cm);
targetLoc = ScreenInfo.xmid + x_loc_pixel;
RNcoordinates = randn(2,1);
dots_targetLoc_coordinates = [targetLoc+(...
    ScreenInfo.px_per_cm.*VSinfo.SD_blob.*RNcoordinates(1,:));...
    ScreenInfo.liftingYaxis+(ScreenInfo.px_per_cm.*...
    VSinfo.SD_yaxis.*RNcoordinates(2,:))];
while 1
    %randomly draw 10 (x,y) coordinates based on the centroid
    RNcoordinates = randn(2,1);
    new_dot_targetLoc_coordinates = [targetLoc+(...
        ScreenInfo.px_per_cm.*VSinfo.SD_blob.*RNcoordinates(1,:));...
        ScreenInfo.liftingYaxis+(ScreenInfo.px_per_cm.*...
        VSinfo.SD_yaxis.*RNcoordinates(2,:))];
    dots_targetLoc_coordinates = [dots_targetLoc_coordinates,new_dot_targetLoc_coordinates];
    %make sure the center of the 10 blobs are aligned with the
    %predetermined location of the test stimulus
    dots_targetLoc_coordinates_shifted = shiftDotClouds(...
        dots_targetLoc_coordinates,x_loc_pixel,ScreenInfo);
    
    %check if they are within the boundaries
    check_withinTheLimit = CheckWithinTheBoundaries(...
        dots_targetLoc_coordinates_shifted,VSinfo.boxSize,ScreenInfo);
    
    %if the generated dots are within boundaries, then pass the
    %coordinates to the function generateDotClouds that gives out the
    %image texture.
    if check_withinTheLimit == 1
        if size(dots_targetLoc_coordinates,2) == VSinfo.num_randomDots
            cmp_stim = generateDotClouds(windowPtr,...
                dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
            break;
        end
    else
        dots_targetLoc_coordinates = dots_targetLoc_coordinates(:,1:end-1);
    end
end

Resp.standard_loc(k,j,i) = ExpInfo.standard_loc(k);
Resp.discrepancy(k,j,i) = Resp.comparison_loc(k,j,i) - Resp.standard_loc(k,j,i);
Resp.staircase(k,j,i) = j;
Resp.order(k,j,i) = ExpInfo.order(k,j,i);

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

if ExpInfo.practice && i == 1; KbWait(-3); end

% blank screen 1
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank1);

if ExpInfo.order(k,j,i) == 1 %1: present the standard first
    
    Screen('DrawTexture',windowPtr,std_stim,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl = Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5)*ExpInfo.IFI);
    
    if ExpInfo.practice && i == 1; KbWait(-3); end
    
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);
    
    Screen('DrawTexture',windowPtr,cmp_stim,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl = Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5)*ExpInfo.IFI);
    
else
    
    Screen('DrawTexture',windowPtr,cmp_stim,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl = Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5)*ExpInfo.IFI);
    
    if ExpInfo.practice && i == 1; KbWait(-3); end
    
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);
    
    Screen('DrawTexture',windowPtr,std_stim,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl = Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr, vbl + (ExpInfo.tStimFrame - 0.5)*ExpInfo.IFI);
    
end

% blank screen 2
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank2);

%% response

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'The first one is more to the right: press 1',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-60,[255 255 255]);
DrawFormattedText(windowPtr, 'The second one is more to the right: press 2',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,[255 255 255]);
Screen('Flip',windowPtr);

tic
while 1
    [~, ~, keyCode] = KbCheck(-1);
    Resp.RT(k,j,i) = toc;
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    elseif keyCode(KbName('1'))
        if ExpInfo.order(k,j,i) == 1 % standard presented first
            Resp.resp(k,j,i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(k,j,i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    elseif keyCode(KbName('2'))
        if ExpInfo.order(k,j,i) == 2 % comparison presented first
            Resp.resp(k,j,i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(k,j,i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    end
end
current_loc = Resp.comparison_loc(k,j,i);
current_resp = Resp.resp(k,j,i);
if rem(j, 2) == 1 && current_resp == -1
    Resp.correct(k,j,i) = 1;
elseif rem(j, 2) == 0 && current_resp == 1
    Resp.correct(k,j,i) = 1;
else
    Resp.correct(k,j,i) = -1;
end

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

%% calcualte location of the comparison stimulus for the next trial

% only update for normal trials
if flag_easy == 0
    
    if i == 1; history = []; else; history = Resp.resp(k,j,1:(i-1)); end
    
    % update location of the comparison stimulus for the next trial
    Resp.comparison_loc(k,j,i+1) = findNextLocation(i,j,current_loc, current_resp, history, ExpInfo);
    
end
end