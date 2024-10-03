
%% Enter experiment info
clear; close all;  rng('Shuffle');

ExpInfo.subjInit = [];
while isempty(ExpInfo.subjInit) == 1
    try ExpInfo.subjInit = input('Participant Initial#: ','s') ;
        ExpInfo.session = input('Session: A/V#: ','s');
        ExpInfo.practice  = input('Main expt: 0; Practice: 1#: ');
    catch
    end
end

 switch ExpInfo.practice
     case 0
         outFileName = sprintf('uniLoc_sub-%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
         ExpInfo.nRep = 30; % number of trial per condition level
         ExpInfo.numBlocks = 4;
     case 1
         outFileName = sprintf('uniLoc_practice_sub-%s_ses-%s', ExpInfo.subjInit, ExpInfo.session);
         ExpInfo.nRep = 5; % number of trial per condition level
         ExpInfo.numBlocks = 2;
 end

 % path control
 curDir = pwd;
[projectDir, ~]  = fileparts(fileparts(curDir));
[git_dir, ~] = fileparts(projectDir);
addpath(genpath(fullfile(git_dir, 'Psychtoolbox-3')))
outDir = fullfile(projectDir, 'data','uniLoc');
if ~exist(outDir,'dir') mkdir(outDir); end

% avoid rewriting data
if exist(fullfile(outDir, [outFileName '.mat']), 'file')
    resp = input('Replace the existing file? Y/N', 's');
    if ~strcmp(resp,'Y')
        disp('Experiment terminated.')
        return
    end
end

if strcmp(ExpInfo.session, 'A')
    if exist('Arduino','var')
        fclose(Arduino);
    end
    Arduino = serial('/dev/cu.usbmodem14301','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
    fopen(Arduino);
    sprintf('Check if system volume is fixed at level 6')
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
ExpInfo.audLevel = [5,7,10,12];
ExpInfo.nLevel = numel(ExpInfo.audLevel);
for tt = 1:ExpInfo.nRep
    ExpInfo.randA(:,tt) = randperm(ExpInfo.nLevel)';
    ExpInfo.randV(:,tt) = randperm(ExpInfo.nLevel)';
end
ExpInfo.randA = reshape(ExpInfo.randA, [], 1)';
ExpInfo.randV = reshape(ExpInfo.randV, [], 1)';

% randomized auditory and visual stimulus location in speaker index
ExpInfo.randAudIdx = ExpInfo.audLevel(ExpInfo.randA);
ExpInfo.randVisIdx = ExpInfo.audLevel(ExpInfo.randV);

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
ExpInfo.nTrials = ExpInfo.nLevel * ExpInfo.nRep;
blocks = linspace(0,ExpInfo.nTrials, ExpInfo.numBlocks+1);
ExpInfo.breakTrials = floor(blocks(2:(end-1)));
ExpInfo.firstTrial = blocks(1:ExpInfo.numBlocks)+1;
ExpInfo.lastTrial = blocks(2:(ExpInfo.numBlocks+1));
ExpInfo.numTrialsPerBlock = ExpInfo.breakTrials(1);

% cost function setup
ExpInfo.maxPoint = 1;
ExpInfo.minPoint = 0.01;
ExpInfo.elbow = ScreenInfo.halfScreenSize*2/4; % in cm
ExpInfo.dropRate = (ExpInfo.maxPoint - ExpInfo.minPoint)/ExpInfo.elbow;

% define durations
ExpInfo.tFixation = 0.3;
ExpInfo.tBlank1 = 0.2;
ExpInfo.tStimFrame = 2;
ExpInfo.tStim = ExpInfo.tStimFrame * ScreenInfo.ifi;
ExpInfo.tIFI = ScreenInfo.ifi;
ExpInfo.tFeedback = 0.8;
ExpInfo.tITI = 0.3;

% dial setup
ExpInfo.dialScaler = 2;
ExpInfo.conf_bar_height = 100;

%% Auditory set up

% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn       = 60;
AudInfo.tStim               = ExpInfo.tStim ; % in sec
duration_gwn                = length(audioSamples)*AudInfo.tStim;
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

% create the bubble visual stimulus
VSinfo.SD_yaxis            = 3; %SD of the blob in cm (vertical)
VSinfo.num_randomDots      = 10; %number of blobs
VSinfo.SD_blob             = 3; %SD of the blob in cm (horizontal)

% draw one blob within the bubbles
VSinfo.width                         = 8; %(pixel) Increasing this value will make the cloud more blurry (arbituary value)
VSinfo.boxSize                       = 15; %This is the box size for each cloud (arbituary value)
VSinfo.maxBrightness                 = 255; %indirectly control contrast
x = 1:1:VSinfo.boxSize; y = x;
[X,Y]                                = meshgrid(x,y);
cloud_temp                           = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
    [VSinfo.width 0; 0 VSinfo.width]);
VSinfo.Cloud                         = reshape(cloud_temp,length(x),length(y)) .* (VSinfo.maxBrightness/max(cloud_temp));

% create background
VSinfo.pblack              = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen          = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture        = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen         = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
VSinfo.blackScreen         = VSinfo.grey_texture;

%% Run the experiment

instruction = ['In the following session, you will be presented \nan auditory or visual stimulus on each trial.','\nAfter the presentation, please use the cursor \nto locate the center of the sound source \n or the center of the bubbles.','\nUse the dial to gain points, and click left to confirm.','\nPress any key to start the unimodal localization task.'];

%start the experiment
c                   = clock;
ExpInfo.start       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5) ,ceil(c(6)));
 
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, instruction ,...
    'center',ScreenInfo.yaxis-300,[255 255 255]);
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

        blockInfo = sprintf('You''ve finished block %i/%i.',idxBlock,ExpInfo.numBlocks);
        scoreInfo1 = sprintf('\nYour cumulative point is %.2f across %i trials.', sum([Resp(1:i).point]), i);
        scoreInfo2 = sprintf('\nYour maximum possible point is %.2f across %i trials.', sum([Resp(1:i).maxPtPossible]), i);
        buttonInfo = '\nPress any button to resume the task.';
        Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        DrawFormattedText(windowPtr, [blockInfo, scoreInfo1, scoreInfo2, buttonInfo],...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
        
    end
end

%% Save sorted data and end the experiment
c  = clock;
ExpInfo.finish  = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
if exist('Arduino','var')
    fclose(Arduino);
end

% sort trials by location level
[~, temp_idx] = sort([Resp(1:end).target_idx]);
sortedResp = Resp(temp_idx);
save(fullfile(outDir,outFileName),'Resp','sortedResp','ExpInfo','ScreenInfo','VSinfo','AudInfo');

%% display leaderoard

if  ExpInfo.practice == 0
    leaderboardText = updateLeaderboardUnimodal(outDir, ExpInfo, Resp);
else
    leaderboardText = 'This is the end of this session. Thank you!';
end
Screen('TextSize',windowPtr,25);
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, leaderboardText,...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
ShowCursor;
Screen('CloseAll');