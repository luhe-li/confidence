function Resp = present_wide_dots(i, j, current_sig, ScreenInfo, ExpInfo, VSinfo, windowPtr, Resp)


%% generate standard stimulus

target_loc = 0; % present at the center
for ff = 1:ExpInfo.n_frame

    while 1

        %randomly draw n_number (x,y) coordinates based on the centroid
        RNcoordinates = randn(2, VSinfo.n_dots);
        dots_targetLoc_coordinates = [target_loc + (...
            ScreenInfo.px_per_cm .* ExpInfo.sig_dot .* RNcoordinates(1,:));...
            ScreenInfo.liftingYaxis + (ScreenInfo.px_per_cm .*...
            VSinfo.SD_yaxis .* RNcoordinates(2,:))];

        %make sure the center of the 10 blobs are aligned with the
        %predetermined location of the test stimulus
        dots_targetLoc_coordinates_shifted = shiftDotClouds(...
            dots_targetLoc_coordinates,target_loc,ScreenInfo);

        %check if they are within the boundaries
        check_withinTheLimit = CheckWithinTheBoundaries(...
            dots_targetLoc_coordinates_shifted,VSinfo.boxSize,ScreenInfo);

        %if the generated dots are within boundaries, then pass the
        %coordinates to the function generateDotClouds that gives out the
        %image texture.
        if check_withinTheLimit == 1
            dotClouds_targetLoc(ff) = generateDotClouds(windowPtr,...
                dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
            break;
        end

    end

end

%% generate comparison stimulus

current_sig = 40;%[0:4:16];
delta_sample_loc = randn(1, ExpInfo.n_frame/2) .* current_sig; % sample target locations per frame
sample_loc = [target_loc + delta_sample_loc, target_loc - delta_sample_loc];
sorted_loc = sort(sample_loc);
sample_loc = sorted_loc([ExpInfo.n_frame/2:-1:1, ExpInfo.n_frame/2+1:ExpInfo.n_frame]);
% r_idx = randperm(ExpInfo.n_frame);
% sample_loc = sample_loc(r_idx);

for ff = 1:ExpInfo.n_frame

    while 1

        %randomly draw n_number (x,y) coordinates based on the centroid
        RNcoordinates = randn(2, VSinfo.n_dots);
        dots_targetLoc_coordinates = [sample_loc(ff) + (...
            ScreenInfo.px_per_cm .* ExpInfo.sig_dot .* RNcoordinates(1,:));...
            ScreenInfo.liftingYaxis + (ScreenInfo.px_per_cm .*...
            VSinfo.SD_yaxis .* RNcoordinates(2,:))];

        %make sure the center of the 10 blobs are aligned with the
        %predetermined location of the test stimulus
        dots_targetLoc_coordinates_shifted = shiftDotClouds(...
            dots_targetLoc_coordinates,sample_loc(ff),ScreenInfo);

        %check if they are within the boundaries
        check_withinTheLimit = CheckWithinTheBoundaries(...
            dots_targetLoc_coordinates_shifted,VSinfo.boxSize,ScreenInfo);

        %if the generated dots are within boundaries, then pass the
        %coordinates to the function generateDotClouds that gives out the
        %image texture.
        if check_withinTheLimit == 1
            cmp_dotClouds_targetLoc(ff) = generateDotClouds(windowPtr,...
                dots_targetLoc_coordinates_shifted,VSinfo,ScreenInfo);
            break;
        end

    end

end


%% trial start

% fixation
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.t_fixation);

% blank screen 1
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.t_blank1);

%1: present standard stimulus first
if ExpInfo.Order(i,j) == 1

    % standard stimulus
    for ff = 1:ExpInfo.n_frame
        for tt = 1: ExpInfo.stay_frame
            Screen('DrawTexture',windowPtr,dotClouds_targetLoc(ff),[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            Screen('Flip',windowPtr);
        end
    end

    % ISI
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(1);

    % comparison stimulus
    for ff = 1:ExpInfo.n_frame
        for tt = 1: ExpInfo.stay_frame
            Screen('DrawTexture',windowPtr,cmp_dotClouds_targetLoc(ff),[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            Screen('Flip',windowPtr);
        end
    end

%2: present comparison stimulus first
else

    % standard stimulus
    for ff = 1:ExpInfo.n_frame
        for tt = 1: ExpInfo.stay_frame
        Screen('DrawTexture',windowPtr,dotClouds_targetLoc(ff),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
        end
    end

    % ISI
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(1);

    % comparison stimulus
    for ff = 1:ExpInfo.n_frame
        for tt = 1: ExpInfo.stay_frame
        Screen('DrawTexture',windowPtr,cmp_dotClouds_targetLoc(ff),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
         Screen('Flip',windowPtr);
        end
    end

end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);







% present stimulus
for ff = 1:ExpInfo.n_frame
Screen('DrawTexture',windowPtr,dotClouds_targetLoc(ff),[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);


DrawFormattedText(windowPtr, 'Debug time...',...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr);




end