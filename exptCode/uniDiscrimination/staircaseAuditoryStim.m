function Resp = staircaseAuditoryStim(i, j, k, ExpInfo, ScreenInfo,AudInfo,Arduino,pahandle,VSinfo,windowPtr,Resp, flag_easy)

% response has the dimension of k_standard_loc x j_staircase x i_trial

if ~exist('flag_easy','Var'); flag_easy = 0; end
if i == ExpInfo.n_trial; flag_easy = 1; end

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

if ExpInfo.order(k,j,i) == 1 %1: present the standard first
    
    % present standard auditory stimulus
    input_on = ['<',num2str(1),':',num2str(ExpInfo.standard_loc),'>'];
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.quietGaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(AudInfo.stimDura);
    input_off = ['<',num2str(0),':',num2str(ExpInfo.standard_loc),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1);
    
    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);
    
    % present comparison auditory stimulus
    input_on = ['<',num2str(1),':',num2str(Resp.comparison_loc(k,j,i)),'>'];
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(AudInfo.stimDura);
    input_off = ['<',num2str(0),':',num2str(Resp.comparison_loc(k,j,i)),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1);

else

    % present comparison auditory stimulus
    input_on = ['<',num2str(1),':',num2str(Resp.comparison_loc(k,j,i)),'>'];
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(AudInfo.stimDura);
    input_off = ['<',num2str(0),':',num2str(Resp.comparison_loc(k,j,i)),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1);

    Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    Screen('Flip',windowPtr);
    WaitSecs(ExpInfo.ISI);
    
    % present standard auditory stimulus
    input_on = ['<',num2str(1),':',num2str(ExpInfo.standard_loc),'>'];
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.quietGaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(AudInfo.stimDura);
    input_off = ['<',num2str(0),':',num2str(ExpInfo.standard_loc),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1);
    
end

% blank screen 2
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank2);

Resp.standard_loc(k,j,i) = mean(ExpInfo.standard_loc(k));
Resp.discrepancy(k,j,i) = Resp.comparison_loc(k,j,i) - Resp.standard_loc(k,j,i);
Resp.staircase(k,j,i) = j;
Resp.order(k,j,i) = ExpInfo.order(k,j,i);

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
    
    if i == 1; history = []; else; history = Resp.resp(k, j,1:(i-1)); end
    
    % update location of the comparison stimulus for the next trial
    Resp.comparison_loc(k,j,i+1) = findNextLocation(i,j,current_loc, current_resp, history, ExpInfo);
    
    % in case of auditory stimulus, make sure that the comparison location
    % is within the speaker range
    Resp.comparison_loc(k,j,i+1) = min(Resp.comparison_loc(k,j,i+1),ExpInfo.n_speaker);
    Resp.comparison_loc(k,j,i+1) = max(Resp.comparison_loc(k,j,i+1),1);
    
end

end