
%% Enter experiment info
clear; close all;  rng('Shuffle');
% 
% ExpInfo.subjInit = [];
% while isempty(ExpInfo.subjInit) == 1
%     try ExpInfo.subjInit = input('Participant Initial#: ','s') ;
%         ExpInfo.session = input('Session: A/V#: ','s');
%         ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
%     catch
%     end
% end

 ExpInfo.subjInit = 'LL';
 ExpInfo.session = 'V';
 ExpInfo.practice  = 1;
        
switch ExpInfo.session
    case 'A'
        ExpInfo.nReliability = 1;
    case 'V'
        ExpInfo.nReliability = 2;
end

switch ExpInfo.practice
    case 1
        outFileName = sprintf('uniLoc_sub%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
        ExpInfo.nRep = 20; % number of trial per condition level
        ExpInfo.numBlocks = 8;
    case 2
        switch ExpInfo.session
            case 'A'
                ExpInfo.nRep = 4;
            case 'V'
                ExpInfo.nRep = 2;          
        end
        outFileName = sprintf('uniLoc_practice_sub%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
        ExpInfo.numBlocks = 2;
end

% path control
curDir = pwd;
[projectDir, ~]  = fileparts(fileparts(curDir));
[git_dir, ~] = fileparts(projectDir);
outDir = fullfile(projectDir, 'data','uniLoc');
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
    % Use INSTRFIND to determine if other instrument objects are connected to the requested device.
    %Arduino = serial('/dev/cu.usbmodemFD131','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
    Arduino = serial('/dev/cu.usbmodem14301','BaudRate',115200);
    fopen(Arduino);
end

%% Screen Setup
PsychDefaultSetup(2);
AssertOpenGL();
GetSecs();
WaitSecs(0.1);
KbCheck();
% ListenChar(2);

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
screens = Screen('Screens');
screenNumber = max(screens);
[windowPtr,rect] = Screen('OpenWindow', screenNumber, [0,0,0]);
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
ScreenInfo.numPixels_perCM = ScreenInfo.xaxis/(ScreenInfo.halfScreenSize*2);

%fixation locations
ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1; 
ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;

%% Experiment set up

% choose auditory locations out of 16 speakers, level/index is speaker
% order (left to right: 1-16)
ExpInfo.audLevel = [5,8,9,12];
ExpInfo.nLevel = numel(ExpInfo.audLevel);
for tt = 1:ExpInfo.nRep
    ExpInfo.randA(:,tt) = randperm(ExpInfo.nLevel)';
    ExpInfo.randV(:,tt) = randperm(ExpInfo.nLevel * ExpInfo.nReliability)';
end
ExpInfo.randA = reshape(ExpInfo.randA, [], 1)';
ExpInfo.randV = reshape(ExpInfo.randV, [], 1)';

% randomized auditory and visual stimulus location in speaker index
ExpInfo.randAudIdx = ExpInfo.audLevel(ExpInfo.randA);
ExpInfo.randVisIdx = ExpInfo.audLevel(ceil(ExpInfo.randV./2));

% location of speakers in CM, visual angle, and pixel
ExpInfo.sittingDistance              = 113.0; %cm
ExpInfo.numSpeaker                   = 16;
ExpInfo.LRmostSpeakers2center        = 65.5; %cm
ExpInfo.LRmostVisualAngle            = (180/pi) * atan(ExpInfo.LRmostSpeakers2center / ...
    ExpInfo.sittingDistance);
ExpInfo.speakerLocCM = linspace(-ExpInfo.LRmostSpeakers2center, ExpInfo.LRmostSpeakers2center, ExpInfo.numSpeaker);
ExpInfo.speakerLocVA = (180/pi) * atan(ExpInfo.speakerLocCM/ExpInfo.sittingDistance);
ExpInfo.speakerLocPixel = round(ExpInfo.speakerLocCM * ScreenInfo.numPixels_perCM);

% auditory locations in different units
ExpInfo.randAudCM     = ExpInfo.speakerLocCM(ExpInfo.randAudIdx);
ExpInfo.randAudVA     = rad2deg(atan(ExpInfo.randAudCM/ExpInfo.sittingDistance));
ExpInfo.randAudPixel  = ExpInfo.randAudCM .* ScreenInfo.numPixels_perCM;

% visual locations in different units
ExpInfo.randVisCM     = ExpInfo.speakerLocCM(ExpInfo.randVisIdx);
ExpInfo.randVisVA     = rad2deg(atan(ExpInfo.randVisCM/ExpInfo.sittingDistance));
ExpInfo.randVisPixel  = ExpInfo.randVisCM .* ScreenInfo.numPixels_perCM;

% split all the trials into blocks
ExpInfo.nTrials = ExpInfo.nLevel * ExpInfo.nRep * ExpInfo.nReliability;
blocks = linspace(0,ExpInfo.nTrials, ExpInfo.numBlocks+1);
ExpInfo.breakTrials = floor(blocks(2:(end-1)));
ExpInfo.firstTrial = blocks(1:ExpInfo.numBlocks)+1;
ExpInfo.lastTrial = blocks(2:(ExpInfo.numBlocks+1));
ExpInfo.numTrialsPerBlock = ExpInfo.breakTrials(1);

% cost function setup
ExpInfo.maxPoint = 100;
ExpInfo.minPoint = 1;
ExpInfo.dropRate = 2;

% define durations
ExpInfo.tFixation = 0.5;
ExpInfo.tBlank1 = 0.3;
ExpInfo.tStimFrame = 9;
ExpInfo.ITI = 0.3;
ExpInfo.tIFI = ScreenInfo.ifi;

%% Auditory set up

% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn       = 100;
AudInfo.stimDura            = ExpInfo.tStimFrame * ExpInfo.tIFI; % in sec
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN       = 1; % too loud for debugging, originally 15
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
                            AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2); % open device

%% audio test / warm-up

if strcmp(ExpInfo.session, 'A')
    for i = 8
    testSpeaker = i;
    input_on = ['<',num2str(1),':',num2str(testSpeaker),'>']; 
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise) 
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(0.1);
    input_off = ['<',num2str(0),':',num2str(testSpeaker),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.5)
    end
end

%% make visual stimuli

% randomize visual reliability by horizontal blob S.D. 
VSinfo.SD_blob(~~rem(ExpInfo.randV,2)) = 8; % in CM
VSinfo.SD_blob(~rem(ExpInfo.randV,2)) = 20; % in CM
VSinfo.numFrames           = ExpInfo.tStimFrame; %for visual stimuli

% create background
VSinfo.pblack              = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen          = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture        = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen         = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);

%% Run the experiment
instruction = ['In the following session, you will be presented \nan auditory or visual stimulus on each trial.','\nAfter the presentation, please use the cursor \nto locate the center of the sound source \n or the center of the visual cloud of dots, not the mode of the cloud.','\nPress key A, S, D, or F to report your confidence in this localization response. ','\nA = Very Low  S = Low  D = High  F = Very High','\nThe key press is used to register localization response.','\nPlease use the whole confidence range.','\nPlease use the same strategy to report your confidence in every session.','\nPress any key to start the unimodal localization task.'];

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

for i = 1:ExpInfo.nTrials
    
    %% present stimuli
    
    if strcmp(ExpInfo.session, 'A')
        SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
        HideCursor;
        Resp(i) = LocalizeAuditoryStim(i, ExpInfo,...
            ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr);
        
    else
        SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
        HideCursor;
          Resp(i)= LocalizeVisualStim(i, ExpInfo,...
            ScreenInfo,VSinfo,windowPtr);
    end
    
    %% save by trial
    save(fullfile(outDir,outFileName),'Resp','ExpInfo','ScreenInfo','VSinfo','AudInfo');
    
    %% add breaks
    if ismember(i,ExpInfo.breakTrials)
        
        Screen('TextSize',windowPtr,30);
        idxBlock = find(ExpInfo.breakTrials==i);
        firstTrial = ExpInfo.firstTrial(idxBlock);
        lastTrial = ExpInfo.lastTrial(idxBlock);

        blockInfo = sprintf('You''ve finished block %i/%i. Please take a break.',idxBlock,ExpInfo.numBlocks);
        Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        DrawFormattedText(windowPtr, blockInfo,...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
        
    end
end

%% Save sorted data and end the experiment
c  = clock;
ExpInfo.finish  = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

if ExpInfo.session == 'V'
    
    [~, temp1] = sort(VSinfo.SD_blob);
    reliSortResp = Resp(temp1);
    reli1resp = reliSortResp(1:ExpInfo.nRep * ExpInfo.nLevel);
    reli2resp = reliSortResp((ExpInfo.nRep * ExpInfo.nLevel+1):end);
    
    [~, temp2] = sort([reli1resp.target_idx]);
    sortedReli1Resp = reli1resp(temp2);
    [~, temp3] = sort([reli2resp.target_idx]);
    sortedReli2Resp = reli2resp(temp3);
    
    save(fullfile(outDir,outFileName),'Resp','reliSortResp','ExpInfo','ScreenInfo','VSinfo','AudInfo','sortedReli1Resp','sortedReli2Resp')
    
else
    % sort trials by location level
    [~, temp2] = sort([Resp(1:end).target_idx]);
    sortedResp = Resp(temp2);
    save(fullfile(outDir,outFileName),'Resp','sortedResp','ExpInfo','ScreenInfo','VSinfo','AudInfo');
    fopen(Arduino);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'End of this session.\nPress any button to exit.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
ShowCursor;
Screen('CloseAll');