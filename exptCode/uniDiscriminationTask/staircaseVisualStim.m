function Resp = staircaseVisualStim(i, j, ExpInfo, ScreenInfo, VSinfo, windowPtr)

%% precompute

% standard stimulus
x_loc_pixel = round(ExpInfo.standard_loc * ExpInfo.px_per_cm);
std_stimTx = generateRippleStim(VSinfo,ScreenInfo,windowPtr, x_loc_pixel);
Resp.standard_loc(i) = ExpInfo.standard_loc;

% comparison stimulus
x_loc_pixel = round(Resp.comparison_loc(i) * ExpInfo.px_per_cm);
cmp_stimTx = generateRippleStim(VSinfo,ScreenInfo,windowPtr, x_loc_pixel);

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
Screen('DrawTexture', windowPtr, std_stimTx(jj),[],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);

%% response

DrawFormattedText(windowPtr, 'The first one is more to the right: press 1',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-60,[255 255 255]);
DrawFormattedText(windowPtr, 'The second one is more to the right: press 2',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,[255 255 255]);
Screen('Flip',windowPtr); WaitSecs(0.1);

tic
while 1
    [~, ~, keyCode] = KbCheck(-1);
    Resp.RT(i) = toc;
    if keyCode(KbName('ESCAPE'))
        sca;
        ShowCursor;
        Screen('CloseAll');
        error('Escape');
    elseif keyCode(KbName('1'))
        if ExpInfo.order(i,j) == 1 % standard presented first
            Resp.resp(i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    elseif keyCode(KbName('2'))
        if ExpInfo.order(i,j) == 2 % standard presented first
            Resp.resp(i) = -1; % participants think the comparison is to the left of the standard
        else
            Resp.resp(i) = 1; % participants think the comparison is to the right of the standard
        end
        break;
    end
end
Screen('Flip',windowPtr);

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

%% calcualte location of the comparison stimulus for the next trial

%1-"LEFT"-up-2-"RIGHT"-down, comparison starts from the LEFT side of the
%standard
current_loc = Resp.comparison_loc(i);
current_resp = Resp.resp(i);
history = Resp.resp(1);

if rem(ExpInfo.condition(i,j),2)==1
    %check whether the participants get the first "LEFT" response yet.
    %If not, the difficulty level will stay the same.
    if i == 1 %when this is the first trial
        %if the subject thinks comparison is on the LEFT side of standard (which is correct)
        if current_resp == -1
            updatedDistance = current_loc + ExpInfo.StepSizes(1); %make it harder
        else
            updatedDistance = current_loc;
        end
    else
        %determine the step size by using the function findPair (find the 
            %number of two consecutive "RIGHT")
          numPairOfRight = findPair([history current_resp],[1 1 -1]);
        
    end
    




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