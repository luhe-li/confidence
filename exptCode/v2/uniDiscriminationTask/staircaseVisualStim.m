function Resp = staircaseVisualStim(i, j, ExpInfo, ScreenInfo, VSinfo, windowPtr, Resp, flag_easy)

if ~exist('flag_easy','Var'); flag_easy = 0; end
if i == ExpInfo.n_trial; flag_easy = 1; end

%% precompute

% standard stimulus
x_loc_pixel = round(ExpInfo.standard_loc * ExpInfo.px_per_cm);
std_stimTx = generateRippleStim(VSinfo,ScreenInfo,windowPtr, x_loc_pixel, 0);
Resp.standard_loc(j,i) = ExpInfo.standard_loc;

% comparison stimulus
x_loc_pixel = round(Resp.comparison_loc(j,i) * ExpInfo.px_per_cm);
cmp_stimTx = generateRippleStim(VSinfo,ScreenInfo,windowPtr, x_loc_pixel, 0);
Resp.discrepancy(j,i) = Resp.comparison_loc(j,i) - Resp.standard_loc(j,i);

Resp.staircase(j,i) = j;
Resp.order(j,i) = ExpInfo.order(j,i);

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

%     for jj = 1:VSinfo.blank_n_frame
%         Screen('DrawTexture', windowPtr, noise_stimTx(jj),[],...
%             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%         Screen('Flip',windowPtr);
%     end

if ExpInfo.order(j,i) == 1 %1: present the standard first

    for jj = 1:VSinfo.numFrames
        Screen('DrawTexture', windowPtr, std_stimTx(jj),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end

    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);

%     for jj = 1:VSinfo.blank_n_frame
%         Screen('DrawTexture', windowPtr, noise_stimTx(jj),[],...
%             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%         Screen('Flip',windowPtr);
%     end

    for jj = 1:VSinfo.numFrames
        Screen('DrawTexture', windowPtr, cmp_stimTx(jj),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end

else

    for jj = 1:VSinfo.numFrames
        Screen('DrawTexture', windowPtr, cmp_stimTx(jj),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end

    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);
%     for jj = 1:VSinfo.blank_n_frame
%         Screen('DrawTexture', windowPtr, noise_stimTx(jj),[],...
%             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%         Screen('Flip',windowPtr);
%     end

    for jj = 1:VSinfo.numFrames
        Screen('DrawTexture', windowPtr, std_stimTx(jj),[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end

end

% blank screen 2
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank2);

%     for jj = 1:VSinfo.blank_n_frame
%         Screen('DrawTexture', windowPtr, noise_stimTx(jj),[],...
%             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%         Screen('Flip',windowPtr);
%     end
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
    Resp.RT(j,i) = toc;
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    elseif keyCode(KbName('1'))
        if ExpInfo.order(j,i) == 1 % standard presented first
            Resp.resp(j,i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(j,i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    elseif keyCode(KbName('2'))
        if ExpInfo.order(j,i) == 2 % comparison presented first
            Resp.resp(j,i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(j,i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    end
end
current_loc = Resp.comparison_loc(j,i);
current_resp = Resp.resp(j,i);
if rem(j, 2) == 1 && current_resp == -1
    Resp.correct(j,i) = 1;
elseif rem(j, 2) == 0 && current_resp == 1
    Resp.correct(j,i) = 1;
else
    Resp.correct(j,i) = -1;
end

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

%% calcualte location of the comparison stimulus for the next trial

% only update for normal trials
if flag_easy == 0
    
    if i == 1; history = []; else; history = Resp.resp(j,1:(i-1)); end
    
    % update location of the comparison stimulus for the next trial
    Resp.comparison_loc(j,i+1) = findNextLocation(i,j,current_loc, current_resp, history, ExpInfo);
    
end
end