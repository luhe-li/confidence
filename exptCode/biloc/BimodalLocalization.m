
% In this task, participants localized visual and auditory stimulic
% presented alone in separate sessions. Each trial started with a fixation
% cross presented straight ahead for 500 ms, followed by 300 ms of blank
% screen. Then, either an auditory or a visual 100 ms-long stimulus was
% presented, followed by 200 ms of blank screen. Next, the response cursor
% appeared, and participants adjusted the horizontal location of the cursor
% to match that of the stimulus. There was no time constraint for the
% response. Visual feedback of the cursor location was provided during its
% adjustment, but error feedback was not provided. After the response, the
% loudspeaker moved to its new location. The inter-trial interval was
% 300 ms.

% Auditory stimuli were presenetd at nine locations, corresponding to
% speaker indices -12:3:12. Each location was tested 20 times
% pseudorandomly, resulting in 180 trials.

%% Enter experiment info
clear; close all;  rng('Shuffle');

ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Participant ID#: ') ;
        ExpInfo.session = input('Session: V1/V2#: ','s');
        ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
    catch
    end
end

switch ExpInfo.practice
    case 1
        outFileName = sprintf('uniLoc_sub%i_ses-%i', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep = 20; % number of trial per condition level
        ExpInfo.numBlocks = 8;
    case 2
        outFileName = sprintf('uniLoc_practice_sub%i_ses-%i', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep = 1; % number of trial per condition level
        ExpInfo.numBlocks = 2;
end

% path control
curDir           = pwd;
[projectDir, ~]  = fileparts(fileparts(curDir));
outDir = fullfile(projectDir, 'data','uniLoc');
if ~exist(outDir,'dir') mkdir(outDir); end
addpath(genpath(PsychtoolboxRoot))

% avoid rewriting data
if exist(fullfile(outDir, [outFileName '.mat']), 'file')
    resp = input('Replace the existing file? Y/N', 's');
    if ~strcmp(resp,'Y')
        disp('Experiment terminated.')
        return
    end
end

% switch between debug mode
ExpInfo.mode  = 1; %input('Experiment mode: 1; Debug mode: 2#: ');
switch ExpInfo.mode
    case 1 % experiment mode
        windowSize = [];
        opacity = 1;
        %Arduino = serial('/dev/cu.usbmodemFD131','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
        Arduino = serial('/dev/cu.usbmodem14101','BaudRate',115200);
        fopen(Arduino);
    case 2 % debug mode
        windowSize = [100 100 1000 600]; % open a smaller window
        opacity = 0.4;
end

%% Screen Setup
PsychDefaultSetup(2);
AssertOpenGL();
GetSecs();
WaitSecs(0.1);
KbCheck();
% ListenChar(2); % silence the keyboard

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
PsychDebugWindowConfiguration([], opacity)
screens = Screen('Screens');
screenNumber = max(screens);
[windowPtr,rect] = Screen('OpenWindow', screenNumber, [0,0,0],windowSize);
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
Screen('TextSize', windowPtr, 30);
Screen('TextFont', windowPtr,'Times');
Screen('TextStyle', windowPtr,1);

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

%% Auditory set up

% open speakers and create sound stimuli
PsychDefaultSetup(2);

% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn       = 10;
AudInfo.stimDura            = 0.1;
AudInfo.tf                  = 400;
AudInfo.intensity           = 0.65;
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, max(timeline_gwn));
AudInfo.intensity_GWN       = 15;
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
    AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
AudInfo.inBetweenGWN        = AudInfo.intensity*AudInfo.GaussianWhiteNoise;
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

%% audio test
PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
PsychPortAudio('Start',pahandle,1,0,0);
WaitSecs(0.1);
input_off = ['<',num2str(0),':',num2str(8),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);

%% make visual stimuli

if strcmp(ExpInfo.session, 'V1')
    VSinfo.SD_blob = 2;
elseif strcmp(ExpInfo.session, 'V2')
    VSinfo.SD_blob = 8;
end
VSinfo.SD_yaxis            = 2; %SD of the blob in cm (vertical)
VSinfo.num_randomDots      = 10; %number of blobs
VSinfo.numFrames           = 6; %for visual stimuli (100 ms)

% create background
VSinfo.pblack                      = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen                  = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture                = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen                 = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);

% draw one blob
VSinfo.width                         = 8; %(pixel) Increasing this value will make the cloud more blurry (arbituary value)
VSinfo.boxSize                       = 11; %This is the box size for each cloud (arbituary value)
x = 1:1:VSinfo.boxSize; y = x;
[X,Y]                                = meshgrid(x,y);
cloud_temp                           = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
    [VSinfo.width 0; 0 VSinfo.width]);
VSinfo.Cloud                         = reshape(cloud_temp,length(x),length(y)) .* (255/max(cloud_temp));

%% Experiment set up

% choose auditory locations out of 16 speakers, in index
ExpInfo.audLevel = [2,4,6,8,9,11,13,15];
ExpInfo.nLevel = numel(ExpInfo.audLevel);
for tt = 1:ExpInfo.nRep
    ExpInfo.randIdx(:,tt) = randperm(ExpInfo.nLevel)';
end
ExpInfo.randIdx = reshape(ExpInfo.randIdx, [], 1)';
ExpInfo.randAudLoc = ExpInfo.audLevel(ExpInfo.randIdx);

% location of speakers in CM, visual angle, and pixel
ExpInfo.sittingDistance              = 113.0; %cm
ExpInfo.numSpeaker                   = 16;
ExpInfo.LRmostSpeakers2center        = 65.5; %cm
ExpInfo.LRmostVisualAngle            = (180/pi) * atan(ExpInfo.LRmostSpeakers2center / ...
                                      ExpInfo.sittingDistance);
ExpInfo.speakerLocCM = linspace(-ExpInfo.LRmostSpeakers2center, ExpInfo.LRmostSpeakers2center, ExpInfo.numSpeaker);
ExpInfo.speakerLocVA = linspace(-ExpInfo.LRmostVisualAngle, ExpInfo.LRmostVisualAngle, ExpInfo.numSpeaker);
ExpInfo.speakerLocPixel = round(ExpInfo.speakerLocCM * ScreenInfo.numPixels_perCM);

% convert target speaker locations in to other units for G.T. visual reference
ExpInfo.loc_idx = ExpInfo.randAudLoc;
ExpInfo.loc_cm  = ExpInfo.speakerLocCM(ExpInfo.randAudLoc);
ExpInfo.loc_deg = rad2deg(atan(ExpInfo.loc_cm/ExpInfo.sittingDistance));
ExpInfo.loc_pixel = ExpInfo.loc_cm .* ScreenInfo.numPixels_perCM;

% split all the trials into blocks
ExpInfo.nTrials = ExpInfo.nLevel * ExpInfo.nRep;
blocks = linspace(0,ExpInfo.nTrials,...
    ExpInfo.numBlocks+1);
ExpInfo.breakTrials       = floor(blocks(2:(end-1)));
ExpInfo.firstTrial     = blocks(1:ExpInfo.numBlocks)+1;
ExpInfo.lastTrial    = blocks(2:(ExpInfo.numBlocks+1));
ExpInfo.numTrialsPerBlock = ExpInfo.breakTrials(1);

% points set up
ExpInfo.maxPoint = 100;
ExpInfo.minPoint = 1; % if enclosed
% maxPoint - droprate * 2 * confidence_radius = minPoint
% we define max 2 * confidence_radius as half of the screen size
ExpInfo.dropRate = (ExpInfo.maxPoint - ExpInfo.minPoint)/ScreenInfo.halfScreenSize;

% define durations
ExpInfo.tFixation = 0.5;
ExpInfo.tBlank1 = 0.3;
ExpInfo.tBlank2 = 0.2;
ExpInfo.tStim = 0.5;
ExpInfo.ITI = 0.3;

%% Run the experiment

%start the experiment
c                   = clock;
ExpInfo.start       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'Press any button to start the unimodal localization task.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
% KbWait(-3);
WaitSecs(1);

for i = 1:ExpInfo.nTrials

    %% present stimuli
    SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
    HideCursor;
    Resp(i) = LocalizeBothStim(i, ExpInfo,...
        ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr);

    %% save by trial
    save(fullfile(outDir,outFileName),'Resp','ExpInfo','ScreenInfo','VSinfo','AudInfo');

    %% add breaks
    if ismember(i,ExpInfo.breakTrials)
        idxBlock = find(ExpInfo.breakTrials==i);
        firstTrial = ExpInfo.firstTrial(idxBlock);
        lastTrial = ExpInfo.lastTrial(idxBlock);
        blockPt = sum([Resp(firstTrial:lastTrial).point]);
        maxPtPossible = lastTrial - firstTrial + 1;
        blockInfo = sprintf('You''ve finished block %i/%i. Please take a break.',idxBlock,ExpInfo.numBlocks);
        pointInfo = sprintf('Your total points of the last block is %.2f (max points possible: %i)',blockPt, maxPtPossible);
        
        Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        DrawFormattedText(windowPtr, blockInfo,...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,...
            [255 255 255]);
        DrawFormattedText(windowPtr, [pointInfo '\nPress any button to resume the experiment.'],...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
    end
end

%% Save sorted data and end the experiment
c  = clock;
ExpInfo.finish  = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

% sort trials by location level
[~, temp] = sort([Resp(1:end).loc_idx]);
sortedResp = Resp(temp);
save(fullfile(outDir,outFileName),'Resp','sortedResp','ExpInfo','ScreenInfo','VSinfo','AudInfo');
 