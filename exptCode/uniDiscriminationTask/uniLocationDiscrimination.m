% Run 2 staircases seperately for A, V1, V2

%% Enter experiment info
clear; close all;  rng('Shuffle');
% 
% ExpInfo.subjInit = [];
% while isempty(ExpInfo.subjInit) == 1
%     try ExpInfo.subjInit = input('Participant Initial#: ','s') ;
%         ExpInfo.session = input('Session: A/V1/V2#: ','s');
%         ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
%     catch
%     end
% end

 ExpInfo.subjInit = 'LL';
 ExpInfo.session = 'V2';
 ExpInfo.practice  = 1;
        
switch ExpInfo.practice
    case 1
        outFileName = sprintf('uniDis_sub%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
    case 2
        outFileName = sprintf('uniDis_practice_sub%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
end


% path control
curDir = pwd;
[projectDir, ~]  = fileparts(fileparts(curDir));
[git_dir, ~] = fileparts(projectDir);
outDir = fullfile(projectDir, 'data','uniDiscrimination');
if ~exist(outDir,'dir') mkdir(outDir); end
addpath(genpath(fullfile(git_dir, 'Psychtoolbox-3')))

% % avoid rewriting data
% if exist(fullfile(outDir, [outFileName '.mat']), 'file')
%     resp = input('Replace the existing file? Y/N', 's');
%     if ~strcmp(resp,'Y')
%         disp('Experiment terminated.')
%         return
%     end
% end

if strcmp(ExpInfo.session, 'A')
    Arduino = serial('/dev/cu.usbmodem14301','BaudRate',115200);
    fopen(Arduino);
    VSinfo.SD_blob = 0; % make sure this variable exists
elseif strcmp(ExpInfo.session, 'V1')
    VSinfo.SD_blob = 4; % cm
elseif strcmp(ExpInfo.session, 'V2')
    VSinfo.SD_blob = 10; % cm
end

%% Screen Setup
PsychDefaultSetup(2);
AssertOpenGL();
GetSecs();
WaitSecs(0.1);
KbCheck();
ListenChar(2);

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
screens = Screen('Screens');
screenNumber = max(screens);
[windowPtr,rect] = Screen('OpenWindow', screenNumber, [0,0,0]); % ,[0,0,800,600]
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
ScreenInfo.screenNumber = screenNumber;
Screen('TextSize', windowPtr, 30);
Screen('TextFont', windowPtr,'Times');
Screen('TextStyle', windowPtr,1);
ScreenInfo.ifi = Screen('GetFlipInterval', windowPtr);

[center(1), center(2)]     = RectCenter(rect);
ScreenInfo.xmid            = center(1); % horizontal center
ScreenInfo.ymid            = center(2); % vertical center
ScreenInfo.backgroundColor = 105;
ScreenInfo.liftingYaxis    = 300;
ScreenInfo.halfScreenSize  = 85; %cm
ScreenInfo.px_per_cm = ScreenInfo.xaxis/(ScreenInfo.halfScreenSize*2);

%fixation locations
ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1; 
ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;

%% Experiment set up

% trial info
ExpInfo.n_trial = 36; % for each staircase
ExpInfo.n_staircase = 2; % one from left and one from right
ExpInfo.n_block = 2;
ExpInfo.n_total_trial = ExpInfo.n_staircase * ExpInfo.n_trial;

% width/distance info
ExpInfo.speaker_cm = 65.5; % cm, left to center
ExpInfo.sitting_dist = 113; % cm
ExpInfo.screen_cm = 170; % cm
ExpInfo.px_per_cm = ScreenInfo.px_per_cm;

% choose speaker/visual location in the main localization experiment
ExpInfo.speaker_idx = 6:11;
ExpInfo.n_speaker = 16;
ExpInfo.speaker_level_cm = linspace(-ExpInfo.speaker_cm, ExpInfo.speaker_cm, ExpInfo.n_speaker);
ExpInfo.sA_cm = ExpInfo.speaker_level_cm(ExpInfo.speaker_idx);
ExpInfo.sV_cm = ExpInfo.sA_cm;

% duration
ExpInfo.n_frame = 10;
ExpInfo.t_fixation = 0.5;
ExpInfo.t_blank1 = 0.2;
ExpInfo.ISI = 0.5;
ExpInfo.ITI = 0.3;
ExpInfo.IFI = ScreenInfo.ifi; 

% split all the trials into blocks
blocks = linspace(0,ExpInfo.n_trial, ExpInfo.n_block+1);
ExpInfo.breakTrials = floor(blocks(2:(end-1)));
ExpInfo.firstTrial = blocks(1:ExpInfo.n_block)+1;
ExpInfo.lastTrial = blocks(2:(ExpInfo.n_block+1));
ExpInfo.numTrialsPerBlock = ExpInfo.breakTrials(1);

% define durations
ExpInfo.tFixation = 0.5;
ExpInfo.tBlank1 = 0.3;
ExpInfo.tStimFrame = 18;
ExpInfo.ISI = 1;
ExpInfo.ITI = 0.3;
ExpInfo.tBlank2 = 0.3;
ExpInfo.tIFI = ScreenInfo.ifi;

% dial
ExpInfo.dialScaler = 2;

% stimulus location for A/V separately
if  strcmp(ExpInfo.session, 'A')
    ExpInfo.standard_loc = [8,9]; % speaker index
    ExpInfo.comparison_loc = 1:ExpInfo.n_speaker;
    ExpInfo.StepSizes = [4,2,1]; % speaker index
else
    ExpInfo.standard_loc = 0; % cm
    ExpInfo.comparison_loc = ExpInfo.sV_cm;
    ExpInfo.StepSizes = [5, 2, 0.5]; % cm
end

%% define staircase conditions

% catch trials are evenly spread across blocks
ExpInfo.n_easy_trial_per_s = 4; % easy trial per staircase
n_trial_w_easy = ExpInfo.n_easy_trial_per_s + ExpInfo.n_trial;

[ExpInfo.condition,ExpInfo.order,...
    Resp.comparison_loc, Resp.standard_loc,...
    Resp.discrepancy, Resp.resp,...
    Resp.RT, Resp.correct] = deal(NaN(ExpInfo.n_staircase,n_trial_w_easy)); 

%------------------------------Conditions----------------------------------
%odd number: starts from leftside of the standard (1-up-2-down) 
%even number: starts from rightside of the standard (2-up-1-down)
for i = 1:n_trial_w_easy
    ExpInfo.condition(:,i) = randperm(ExpInfo.n_staircase,ExpInfo.n_staircase);
end 

Resp.comparison_loc(1,1) = min(ExpInfo.comparison_loc);
Resp.comparison_loc(2,1) = max(ExpInfo.comparison_loc);

%---------------------------------ExpInfo.order------------------------------------
%1: present the standard first
%2: present the comparison first
order(ExpInfo.condition==1) = reshape(Shuffle(repmat([1,2],[1, n_trial_w_easy/2])),[2, n_trial_w_easy/2]);
order(ExpInfo.condition==2) = reshape(Shuffle(repmat([1,2],[1, n_trial_w_easy/2])),[2, n_trial_w_easy/2]);
order = reshape(order, [2, n_trial_w_easy]); % row 1: trials of staircase1; row 2: trials of staircase2

%-------------------------------Response-----------------------------------
%-1: participants think the comparison is to the left of the standard
%1: participants think the comparison is to the right of the standard

%% Add easy trials to calculate lapse rate

trial_slc  = [];
for i = 1:ExpInfo.n_staircase
    trialIdx_staircase = reshape(i:ExpInfo.n_staircase:ExpInfo.n_total_trial, ...
        [ExpInfo.n_trial/ExpInfo.n_block, ExpInfo.n_block]);
    for j = 1:ExpInfo.n_block
        trial_slc = [trial_slc; trialIdx_staircase(randperm(ExpInfo.n_trial/ExpInfo.n_block,...
            ExpInfo.n_easy_trial_per_s/ExpInfo.n_block), j)];
    end
end
ExpInfo.easy_trial = reshape(trial_slc,[ExpInfo.n_easy_trial_per_s, ExpInfo.n_staircase])';
% %To see whether the catch trials are evenly spread out. Do:
% T = zeros(ExpInfo.n_staircase, ExpInfo.n_trial);
% T(trial_slc) = 1; imagesc(T);

% set parameters for easy trials as the last few trials
ExpInfo.easy_idx = ExpInfo.n_trial+1:ExpInfo.n_trial+ExpInfo.n_easy_trial_per_s;
Resp.comparison_loc(1, ExpInfo.easy_idx) = min(ExpInfo.comparison_loc);
Resp.comparison_loc(2, ExpInfo.easy_idx) = max(ExpInfo.comparison_loc);

%% Auditory set up

% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
AudInfo.stimDura            = ExpInfo.tStimFrame * ExpInfo.IFI; % in sec
standardFrequency_gwn       = 1/(ExpInfo.tStimFrame * ExpInfo.IFI);%100;
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN       = 1;
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
    AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
AudInfo.quietGaussianWhiteNoise  =0.5*[AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
    AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn]; % half the amplitude for two speakers of the standard stimulus
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2); % open device

%% audio test / warm-up

if strcmp(ExpInfo.session, 'A')
    for i = 8
    testSpeaker = i;
    input_on = ['<',num2str(1),':',num2str(testSpeaker),'>']; 
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise) ;
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(0.5);
    input_off = ['<',num2str(0),':',num2str(testSpeaker),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1)
    end
end

%% make visual stimuli

% define ripple stimulus 
VSinfo.height = 200; % pixel
VSinfo.noise_sd = 30; % pixel
VSinfo.stim_sd = VSinfo.SD_blob.* ExpInfo.px_per_cm; % doesnt change by trial
VSinfo.wBack = 0.6; % [0, 1]
VSinfo.pContrast = 0.6; % [0, 1]

% create background
VSinfo.numFrames           = ExpInfo.tStimFrame; %for visual stimuli
VSinfo.pblack              = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen          = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture        = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen         = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);

%% Run the experiment
instruction = ['\nPress any key to start the unimodal location discrimination task.'];

%start the experiment
c                   = clock;
ExpInfo.start       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5) ,ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, instruction ,...
    'center',ScreenInfo.yaxis-500,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
WaitSecs(1);

for i = 1:ExpInfo.n_trial

    for ss = 1:ExpInfo.n_staircase
        %% present stimulus

        % order by condition index
        j = ExpInfo.condition(ss,i);
        ExpInfo.order(j,i) = order(ss,i);
        
        % stimulus
        SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
        HideCursor;
        if strcmp(ExpInfo.session, 'A')
            Resp = staircaseAuditoryStim(i, j, ExpInfo,...
                ScreenInfo,AudInfo,Arduino,pahandle,VSinfo,windowPtr,Resp);
        else
            Resp = staircaseVisualStim(i, j, ExpInfo,...
                ScreenInfo,VSinfo,windowPtr,Resp);
        end
        
        %% inserted an easy trial for calculating the lapse rate later
        
        i_total = (ss-1)*ExpInfo.n_trial + i;
        if ismember(i_total, ExpInfo.easy_trial)
            
            flag_easy = 1;
            [j_easy, temp_i_easy] = find(ExpInfo.easy_trial==i_total);
            i_easy = ExpInfo.easy_idx(temp_i_easy);
            
            if strcmp(ExpInfo.session, 'A')
                Resp = staircaseAuditoryStim(i_easy, j_easy, ExpInfo,...
                    ScreenInfo,AudInfo,Arduino,pahandle,VSinfo,windowPtr,Resp,flag_easy);
            else
                Resp = staircaseVisualStim(i_easy, j_easy, ExpInfo,...
                    ScreenInfo,VSinfo,windowPtr, Resp, flag_easy);
            end
            
        end
        
        %% save by trial
        
        save(fullfile(outDir,outFileName),'Resp','ExpInfo','ScreenInfo','VSinfo','AudInfo');
        
        %% add breaks
        if ismember(i,ExpInfo.breakTrials) && ss == 2
            
            Screen('TextSize',windowPtr,30);
            idxBlock = find(ExpInfo.breakTrials==i);
            
            blockInfo = sprintf('You''ve finished block %i/%i. Please take a break.',idxBlock, ExpInfo.n_block);
            Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            DrawFormattedText(windowPtr, blockInfo,...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
                [255 255 255]);
            Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
            
        end
        
    end
    
end
%% Save sorted data and end the experiment
c  = clock;
ExpInfo.finish  = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
% 
% if ExpInfo.session == 'V'
%     
%     [~, temp1] = sort(VSinfo.SD_blob);
%     reliSortResp = Resp(temp1);
%     reli1resp = reliSortResp(1:ExpInfo.nRep * ExpInfo.nLevel);
%     reli2resp = reliSortResp((ExpInfo.nRep * ExpInfo.nLevel+1):end);
%     
%     [~, temp2] = sort([reli1resp.target_idx]);
%     sortedReli1Resp = reli1resp(temp2);
%     [~, temp3] = sort([reli2resp.target_idx]);
%     sortedReli2Resp = reli2resp(temp3);
%     
%     save(fullfile(outDir,outFileName),'Resp','reliSortResp','ExpInfo','ScreenInfo','VSinfo','AudInfo','sortedReli1Resp','sortedReli2Resp')
%     
% else
%     % sort trials by location level
%     [~, temp2] = sort([Resp(1:end).target_idx]);
%     sortedResp = Resp(temp2);
%     save(fullfile(outDir,outFileName),'Resp','sortedResp','ExpInfo','ScreenInfo','VSinfo','AudInfo');
%     fopen(Arduino);
% end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'End of this session.\nPress any button to exit.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);

KbWait(-3);
ShowCursor;
Screen('CloseAll');