
%% Enter experiment info
clear; close all;  rng('Shuffle');

ExpInfo.subjID                       = [];
while isempty(ExpInfo.subjID)        == 1
    try ExpInfo.subjID                   = input('Participant ID#: ') ;
        ExpInfo.session                      = input('Session#: ');
        ExpInfo.practice                     = input('Main expt: 1; Practice: 2#: ');
    catch
    end
end

switch ExpInfo.practice
    case 1
        outFileName                          = sprintf('biLoc_sub%i_ses%i', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep                         = 10; % number of trial per condition level
        ExpInfo.numBlocks                    = 8;
    case 2
        outFileName                          = sprintf('biLoc_practice_sub%i_ses%i', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep                         = 2; % number of trial per condition level
        ExpInfo.numBlocks                    = 2;
end

% path control
curDir                               = pwd;
[projectDir, ~]                      = fileparts(fileparts(curDir));
outDir                               = fullfile(projectDir, 'data','biLoc');
if ~exist(outDir,'dir') mkdir(outDir); end
addpath(genpath(PsychtoolboxRoot))

% avoid rewriting data
if exist(fullfile(outDir, [outFileName '.mat']), 'file')
    resp                                 = input('Replace the existing file? Y/N', 's');
    if ~strcmp(resp,'Y')
        disp('Experiment terminated.')
        return
    end
end

% switch between debug mode
ExpInfo.mode                         = 1; %input('Experiment mode: 1; Debug mode: 2#: ');
switch ExpInfo.mode
    case 1 % experiment mode
        windowSize                           = [];
        opacity                              = 1;
        % make sure this value matches with the baudrate in the arduino code
        Arduino                              = serial('/dev/cu.usbmodem14101','BaudRate',115200);
        fopen(Arduino);
    case 2 % debug mode
        windowSize                           = [100 100 1000 600]; % open a smaller window
        opacity                              = 0.4;
end

%% Screen Setup
PsychDefaultSetup(2);
AssertOpenGL();
GetSecs();
WaitSecs(0.1);
KbCheck();
ListenChar(1); % change from 2 to 1 because if 2 then NO keyboard will be usable after quitting with escape.
% STOP CHANGING THIS TO TWO 
% S T O P!    D O I N Gsssssssssss!    T H A T!
% ListenChar(        ONE!        );
% There isn't even a keyboard with a index of two what are you even listening to 
% DO NOT CHANGE THIS TO TWO

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
PsychDebugWindowConfiguration([], opacity)
screens                              = Screen('Screens');
screenNumber                         = max(screens);
[windowPtr,rect]                     = Screen('OpenWindow', screenNumber, [0,0,0],windowSize);
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
Screen('TextSize', windowPtr, 30);
Screen('TextFont', windowPtr,'Times');
Screen('TextStyle', windowPtr,1);
ScreenInfo.ifi                       = Screen('GetFlipInterval', windowPtr);

[center(1), center(2)]               = RectCenter(rect);
ScreenInfo.xmid                      = center(1); % horizontal center
ScreenInfo.ymid                      = center(2); % vertical center
ScreenInfo.backgroundColor           = 105;
ScreenInfo.liftingYaxis              = 300;
ScreenInfo.halfScreenSize            = 85; %cm
ScreenInfo.numPixels_perCM           = ScreenInfo.xaxis/(ScreenInfo.halfScreenSize*2);

%fixation locations
ScreenInfo.x1_lb                     = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
ScreenInfo.x1_ub                     = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
ScreenInfo.y1_lb                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
ScreenInfo.y1_ub                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
ScreenInfo.y2_lb                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
ScreenInfo.y2_ub                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;

%% Experiment set up

% choose audiovisual locations out of 16 speakers, in index
ExpInfo.audIdx                       = [6,8,9,11];
ExpInfo.visIdx                       = [6,8,9,11];
ExpInfo.cueIdx                       = [1,2]; % 1 = A, 2 = V
ExpInfo.visReliIdx                   = [8,20];%[10,28]; % std in centimeters
ExpInfo.cue                          = {'A','V'};
ExpInfo.avIdx                        = combvec(ExpInfo.audIdx, ExpInfo.visIdx, ExpInfo.cueIdx, ExpInfo.visReliIdx);
ExpInfo.nLevel                       = size(ExpInfo.avIdx, 2);
for tt                               = 1:ExpInfo.nRep
    ExpInfo.randIdx(:,tt)                = randperm(ExpInfo.nLevel)';
end
ExpInfo.randIdx                      = reshape(ExpInfo.randIdx, [], 1)';
ExpInfo.randAVIdx                    = ExpInfo.avIdx(:,ExpInfo.randIdx);
VSinfo.SD_blob                       = ExpInfo.randAVIdx(4,:);

% location of speakers in CM, visual angle, and pixel
ExpInfo.sittingDistance              = 113; %cm
ExpInfo.numSpeaker                   = 16;
ExpInfo.LRmostSpeakers2center        = 65.5; %cm
ExpInfo.LRmostVisualAngle            = (180/pi) * atan(ExpInfo.LRmostSpeakers2center / ...
                                      ExpInfo.sittingDistance);
ExpInfo.speakerLocCM                 = linspace(-ExpInfo.LRmostSpeakers2center, ExpInfo.LRmostSpeakers2center, ExpInfo.numSpeaker);
ExpInfo.speakerLocVA                 = linspace(-ExpInfo.LRmostVisualAngle, ExpInfo.LRmostVisualAngle, ExpInfo.numSpeaker);
ExpInfo.speakerLocPixel              = round(ExpInfo.speakerLocCM * ScreenInfo.numPixels_perCM);

% auditory locations in different units
ExpInfo.randAudIdx    = ExpInfo.randAVIdx(1,:);
ExpInfo.randAudCM     = ExpInfo.speakerLocCM(ExpInfo.randAudIdx);
ExpInfo.randAudVA     = rad2deg(atan(ExpInfo.randAudCM/ExpInfo.sittingDistance));
ExpInfo.randAudPixel  = ExpInfo.randAudCM .* ScreenInfo.numPixels_perCM;

% visual locations in different units
ExpInfo.randVisIdx    = ExpInfo.randAVIdx(2,:);
ExpInfo.randVisCM     = ExpInfo.speakerLocCM(ExpInfo.randVisIdx);
ExpInfo.randVisVA    = rad2deg(atan(ExpInfo.randVisCM/ExpInfo.sittingDistance));
ExpInfo.randVisPixel = ExpInfo.randVisCM .* ScreenInfo.numPixels_perCM;

% convert visual locations from index to perceptually matching pixel
load([sprintf('AVbias_sub%i', ExpInfo.subjID) '.mat'])

x = ExpInfo.speakerLocPixel(ExpInfo.randAudIdx);
coefsA = squeeze(Transfer.PxCoeff(1, :));
fitRA = x .* coefsA(2) + coefsA(1);
fitSV = fitRA - ScreenInfo.xmid; 
ExpInfo.targetPixel                  = unique(fitSV);
[~, ~, ic]                           = unique(ExpInfo.randAVIdx(2,:));
ExpInfo.randVisPixel                 = ExpInfo.targetPixel(ic');

% split all the trials into blocks
if ExpInfo.practice == 1
    ExpInfo.nTrials                      = ExpInfo.nLevel * ExpInfo.nRep;
elseif ExpInfo.practice == 2
    ExpInfo.nTrials                      = 20;
end
blocks                               = linspace(0,ExpInfo.nTrials,...
    ExpInfo.numBlocks+1);
ExpInfo.breakTrials                  = floor(blocks(2:(end-1)));
ExpInfo.firstTrial                   = blocks(1:ExpInfo.numBlocks)+1;
ExpInfo.lastTrial                    = blocks(2:(ExpInfo.numBlocks+1));
ExpInfo.numTrialsPerBlock            = ExpInfo.breakTrials(1);

% define durations
ExpInfo.tFixation                    = 0.5;
ExpInfo.tBlank1                      = 0.3;
ExpInfo.tStimFrame = 3; % in frame
ExpInfo.ITI = 0.3;
ExpInfo.tIFI = ScreenInfo.ifi;

% ExpInfo.frameStim                    = round(AudInfo.stimDura * 60);
% ExpInfo.tStim                        = ExpInfo.frameStim * (1/60);
% ExpInfo.ITI                          = 0.3;

%% Auditory set up

% get correct sound card
InitializePsychSound
devices                              = PsychPortAudio('GetDevices');
our_device                           = devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                           = 44100;
audioSamples                         = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn                = 100;
AudInfo.stimDura                     = ExpInfo.tStimFrame * ExpInfo.tIFI; % in sec
duration_gwn                         = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                         = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn                       = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn                     = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN                = 1; % too loud for debugging, orginally 15
AudInfo.GaussianWhiteNoise           = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
    AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
pahandle                             = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

%% audio test / warm-up
duration_warm                = length(audioSamples);
timeline_warm                = linspace(1,duration_warm,duration_warm);
sineWindow_warm              = sin(standardFrequency_gwn/2*2*pi*timeline_warm/AudInfo.fs);
carrierSound_warm            = randn(1, numel(timeline_warm));

AudInfo.WarmupWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_warm.*carrierSound_warm; AudInfo.intensity_GWN.*sineWindow_warm.*carrierSound_warm];


testSpeaker = 8;
input_on = ['<',num2str(1),':',num2str(testSpeaker),'>']; %arduino takes input in this format
fprintf(Arduino,input_on);
PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
PsychPortAudio('Start',pahandle,1,0,0);
WaitSecs(1);
input_off = ['<',num2str(0),':',num2str(testSpeaker),'>'];
fprintf(Arduino,input_off);
PsychPortAudio('Stop',pahandle);
WaitSecs(1);

%% make visual stimuli

VSinfo.SD_yaxis                      = 5; %SD of the blob in cm (vertical)
VSinfo.num_randomDots                = 10; %number of blobs
VSinfo.numFrames                     = 3; %for visual stimuli
VSinfo.numFramesMasker               = 30; %for mask

% create background
VSinfo.pblack                        = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen                    = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture                  = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen                   = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);

% white noise background
VSinfo.GWNnumPixel                   = 4; % 4 pixels will have the same color
VSinfo.GWNnumFrames                  = 10; 
VSinfo.gwn_texture                   = generateNoisyBackground(VSinfo,ScreenInfo,windowPtr);

% draw one blob
VSinfo.width                         = 8; %(pixel) Increasing this value will make the cloud more blurry (arbituary value)
VSinfo.boxSize                       = 15; %This is the box size for each cloud (arbituary value)
VSinfo.maxBrightness                 = 128; %indirectly control contrast
x                                    = 1:1:VSinfo.boxSize; y = x;
[X,Y]                                = meshgrid(x,y);
cloud_temp                           = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
    [VSinfo.width 0; 0 VSinfo.width]);
VSinfo.Cloud                         = reshape(cloud_temp,length(x),length(y)) .* (VSinfo.maxBrightness/max(cloud_temp));

%% Run the experiment

%start the experiment
c                                    = clock;
ExpInfo.start                        = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'Press any button to start the bimodal localization task.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
WaitSecs(1);

for i                                = 1:ExpInfo.nTrials

    %% present stimuli
    SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
    HideCursor;
    Resp(i)                              = LocalizeBothStim(i, ExpInfo,...
        ScreenInfo,AudInfo,VSinfo,Arduino,pahandle,windowPtr);

    %% save by trial
    save(fullfile(outDir,outFileName),'Resp','ExpInfo','ScreenInfo','VSinfo','AudInfo');

    %% add breaks
    if ismember(i,ExpInfo.breakTrials)
        idxBlock                             = find(ExpInfo.breakTrials==i);
        firstTrial                           = ExpInfo.firstTrial(idxBlock);
        lastTrial                            = ExpInfo.lastTrial(idxBlock);

        blockInfo = sprintf('You''ve finished block %i/%i. Please take a break.',idxBlock,ExpInfo.numBlocks);
     
        Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        DrawFormattedText(windowPtr, blockInfo,...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,...
            [255 255 255]);
        DrawFormattedText(windowPtr, '\nPress any button to resume the experiment.',...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
    end
end

%% Save sorted data and end the experimentz
c                                    = clock;
ExpInfo.finish                       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, 'End of this session.\nPress any button to exit.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
ShowCursor;
Screen('CloseAll');