
% In this task, participants localized visual and auditory stimuli
% presented alone in separate sessions. Each trial started with a fixation
% cross presented straight ahead for 500 ms, followed by 700 ms of blank
% screen. Then, either an auditory or a visual 100 ms-long stimulus was
% presented, followed by 700 ms of blank screen. Next, the response cursor
% appeared, and participants adjusted the horizontal location of the cursor
% to match that of the stimulus. There was no time constraint for the
% response. Visual feedback of the cursor location was provided during its
% adjustment, but error feedback was not provided. After the response, the
% loudspeaker moved to its new location. The inter-trial interval was
% jittered around 1000 ms.

% Auditory stimuli were presenetd at nine locations, corresponding to
% speaker indices -12:3:12. Each location was tested 20 times
% pseudorandomly, resulting in 180 trials.

%% Enter experiment info
clear; close all; clc; rng('Shuffle');

ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Please enter participant ID#: ') ; %'s'
        ExpInfo.session = input('Please enter session#: ');
        ExpInfo.practice  = input('Experiment mode: 1; Practice mode: 2#: ');
    catch
    end
end

switch ExpInfo.practice
    case 1 
        outFileName = sprintf('uniçLoc_sub%s_ses%s', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep = 20; % number of trial per condition level
    case 2 
        outFileName = sprintf('uniLoc_practice_sub%s_ses%s', ExpInfo.subjID, ExpInfo.session);
        ExpInfo.nRep = 1; % number of trial per condition level
end

% path control
outDir = fullfile(pwd, 'data');
if ~exist(outDir,'dir') mkdir(outDir); end
addpath(genpath('/e/3.3/p3/hong/Desktop/Project5/Psychtoolbox'));
addpath(genpath(PsychtoolboxRoot))

% avoid rewriting data
if exist(fullfile(outDir, [outFileName '.mat']), 'file')
    resp                                 = input('To replace the existing file, press y', 's');
    if ~strcmp(resp,'y')
        disp('Experiment terminated.')
        return
    end
end

% switch between debug mode
ExpInfo.mode  = input('Experiment mode: 1; Debug mode: 2#: ');
switch ExpInfo.mode
    case 1 % experiment mode
        windowSize = [];
        opacity = 1;
        HideCursor();
        Arduino = serial('/dev/cu.usbmodemFD131','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
        fopen(Arduino);
    case 2 % debug mode
        windowSize = [100 100 1000 600]; % open a smaller window
        opacity = 0.4;
end

%% Experiment set up

% choose auditory locations out of 31 speakers
audLevel = 4:3:28;
nLevel = numel(audLevel);
audLoc = repmat(audLevel, [1,ExpInfo.nRep]);
randIdx = randperm(numel(audLoc));
randAudLoc = audLoc(randIdx);


%% Screen Setup
PsychDefaultSetup(2);
AssertOpenGL();
GetSecs();
WaitSecs(0.1);
KbCheck();
ListenChar(2); % silence the keyboard

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
PsychDebugWindowConfiguration([], opacity)
screens = Screen('Screens');
screenNumber = max(screens);
[windowPtr,rect] = Screen('OpenWindow', screenNumber, [0,0,0],windowSize);
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
Screen('TextSize', windowPtr, 35) ;   
Screen('TextFont', windowPtr,'Times');
Screen('TextStyle', windowPtr,1); 

[center(1), center(2)]     = RectCenter(rect);
ScreenInfo.xmid            = center(1); % horizontal center
ScreenInfo.ymid            = center(2); % vertical center
ScreenInfo.backgroundColor = 105;
ScreenInfo.numPixels_perCM = 7.5;
ScreenInfo.liftingYaxis    = 300; 

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
our_device=devices(3).DeviceIndex;

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

%% Specify the locations of the auditory stimulus
% store the differences between standard and test four locations
AudInfo.Distance        = round(ExpInfo.matchedAloc,2); %in deg
AudInfo.InitialLoc      = -7.5;
AudInfo.numLocs         = length(AudInfo.Distance); 
AudInfo.numTrialsPerLoc = 30; %1 for practice, 30 for the real experiment
AudInfo.numTotalTrials  = AudInfo.numTrialsPerLoc * AudInfo.numLocs;
%the order of presenting V/A (if 1: A; %if 2: V)
%Even when the trial is V, we still move the speaker
ExpInfo.order_VSnAS     = [];
for i = 1:AudInfo.numTotalTrials
    ExpInfo.order_VSnAS = [ExpInfo.order_VSnAS, randperm(2,2)]; 
end
%shuffle auditory locations
rand_indices = [];
for i = 1:AudInfo.numTrialsPerLoc; rand_indices = [rand_indices, randperm(4,4)]; end

AudInfo.trialConditions = NaN(1,AudInfo.numTotalTrials*2);
AudInfo.trialConditions(ExpInfo.order_VSnAS==1) = AudInfo.Distance(rand_indices); 
AudInfo.trialConditions(ExpInfo.order_VSnAS==2) = Shuffle(AudInfo.Distance(rand_indices)); 

%1st row: location (in deg)
%2nd row: responses (in deg)
%3rd row: response time
AudInfo.data      = zeros(3,AudInfo.numTotalTrials);
AudInfo.data(1,:) = AudInfo.trialConditions(ExpInfo.order_VSnAS==1); %in deg
%create a matrix whose 1st column is initial position and 2nd column is final position
AudInfo.Location_Matrix = [[AudInfo.InitialLoc, AudInfo.trialConditions(1:end-1)]',...
                                 AudInfo.trialConditions'];
[AudInfo.locations_wMidStep,AudInfo.moving_locations_steps, AudInfo.totalSteps] =...
    randomPath2(AudInfo.Location_Matrix, ExpInfo); %1 midpoint, 2 steps
AudInfo.waitTime         = 0.5; 

%% define the locatons of visual stimuli
VSinfo.initialDistance  = -12:8:12; %in deg
VSinfo.numLocs          = length(VSinfo.initialDistance);
VSinfo.numTrialsPerLoc  = AudInfo.numTrialsPerLoc;
VSinfo.numFrames        = 6;
VSinfo.width            = 401; %(pixel) Increasing this value will make the cloud more blurry
VSinfo.boxSize          = 201; %This is the box size for each cloud.
VSinfo.intensity        = 10; %This determines the height of the clouds. Lowering this value will make
                                %them have lower contrast
%set the parameters for the visual stimuli, which consist of 10 gaussian blobs
VSinfo.blackScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
x                  = 1:1:VSinfo.boxSize; y = x;
[X,Y]              = meshgrid(x,y);
cloud              = 1e2.*mvnpdf([X(:) Y(:)],[median(x) median(y)],...
                        [VSinfo.width 0; 0 VSinfo.width]);
VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(cloud,length(x),length(y)); 

%shuffle visual locations
rand_indices = [];
for i = 1:VSinfo.numTrialsPerLoc; rand_indices = [rand_indices, randperm(4,4)]; end
VSinfo.trialConditions = VSinfo.initialDistance(rand_indices); 

%date_easyTrials stores all the data for randomly inserted easy trials
%1st row: the target will appear at either of the four locations (in deg)
%2nd row: the target location (in cm)
%3rd row: response (in deg)
%4th row: Response time
VSinfo.numTotalTrials = VSinfo.numTrialsPerLoc*VSinfo.numLocs;
VSinfo.data           = zeros(4,VSinfo.numTotalTrials);
VSinfo.data(1,:)      = VSinfo.trialConditions;
VSinfo.data(2,:)      = tan(deg2rad(VSinfo.data(1,:))).*ExpInfo.sittingDistance;

%% specify the experiment informations
ExpInfo.numBlocks         = 4;
blocks                    = linspace(0,AudInfo.numTotalTrials+VSinfo.numTotalTrials,...
                            ExpInfo.numBlocks+1); 
%split all the trials into 4 blocks
ExpInfo.breakTrials       = floor(blocks(2:(end-1)));
ExpInfo.numTrialsPerBlock = ExpInfo.breakTrials(1);

%% Run the experiment by calling the function InterleavedStaircase
%start the experiment
DrawFormattedText(windowPtr, 'Press any button to start the unimodal localization task.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3); WaitSecs(1);
Screen('Flip',windowPtr);

for i = 1:(AudInfo.numTotalTrials+VSinfo.numTotalTrials) 
    if ExpInfo.order_VSnAS(i) == 1 %auditory trial
        jj = sum(ExpInfo.order_VSnAS(1:i)==1);
        [AudInfo.data(2,jj), AudInfo.data(3,jj)] = PresentAuditoryStimulus(i,...
            ExpInfo,ScreenInfo,AudInfo,motorArduino,noOfSteps,pahandle,windowPtr);
        
    else %visual trial
        ii = sum(ExpInfo.order_VSnAS(1:i)==2);
        [VSinfo.data(3,ii),VSinfo.data(4,ii)] = PresentVisualStimulus(i,...
            ii,ExpInfo,ScreenInfo,VSinfo,AudInfo,motorArduino,noOfSteps,...
            pahandle,windowPtr);
    end
    
    %add breaks     
    if ismember(i,ExpInfo.breakTrials)
        DrawFormattedText(windowPtr, 'You''ve finished one block. Please take a break.',...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis-30,...
            [255 255 255]);
        DrawFormattedText(windowPtr, 'Press any button to resume the experiment.',...
            'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,...
            [255 255 255]);
        Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
        Screen('Flip',windowPtr); WaitSecs(0.5);
    end  
end

%% Move the arduino to the leftmost place
if abs(AudInfo.InitialLoc - AudInfo.Location_Matrix(end))>=0.01
    steps_goingBack = round((tan(deg2rad(AudInfo.InitialLoc))-...
                      tan(deg2rad(AudInfo.Location_Matrix(end))))*...
                      ExpInfo.sittingDistance/3,2);
    if steps_goingBack < 0 
        fprintf(motorArduino,['%c','%d'], ['p',noOfSteps*abs(steps_goingBack)]);
    else
        fprintf(motorArduino,['%c','%d'], ['n',noOfSteps*abs(steps_goingBack)]);
    end 
end
delete(motorArduino)

%% Save data and end the experiment
Unimodal_localization_data = {ExpInfo,ScreenInfo,VSinfo,AudInfo};
save(outFileName,'Unimodal_localization_data');
Screen('CloseAll');
