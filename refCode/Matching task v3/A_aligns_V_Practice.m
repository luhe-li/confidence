%In this 2-IFC task, participants judged the location of an auditory 
%stimulus relative to that of a visual standard stimulus. The two stimuli 
%were presented sequentially, in pseudo-randomized order. Each trial 
%started with the presentation of a fixation cross located at the center of
%the screen for 500 ms, followed by a blank screen lasting for 1000 ms. 
%Then, the first, 100 ms-long stimulus was presented, followed by a blank 
%screen presented for 1000 ms. After the first stimulus presentation, the 
%same events (fixation, blank screen, second stimulus) were repeated. At 
%the end of each trial, a response probe was displayed. Participants 
%indicated by button press whether the auditory stimulus was located to 
%the left or right of the visual stimulus. Feedback was not provided 
%(Fig. 2A). The inter-trial interval was approximately 3000 ms. The long 
%period was needed to move the loudspeaker to its new position for the next 
%trial.

%This is just a practice

%% Enter subject's name
clear all; close all; clc
rng('shuffle');
%enter subject's ID
ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Please enter participant ID#: ') ; %'s'
    catch
    end
end
%Create a file name that stores all the data of the visual JND experiment
out1FileName = ['A_aligns_V_practice_sub' num2str(ExpInfo.subjID)];

%% Screen Setup 
%Canvas size = 53.5" x 40"= 135.8cm x 101.6cm
%Screen size by the project = 1024 pixels x 768 pixels
%so each centimeter has 7.54 pixels horizontally and 7.56 vertically
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
[windowPtr,rect] = Screen('OpenWindow', 0, [0,0,0]);
%[windowPtr,rect] = Screen('OpenWindow', 0, [0,0,0],[100 100 1000 680]); % for testing
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
Screen('TextSize', windowPtr, 35) ;   
Screen('TextFont',windowPtr,'Times');
Screen('TextStyle',windowPtr,1); 

[center(1), center(2)]      = RectCenter(rect);
ScreenInfo.xmid             = center(1); % horizontal center
ScreenInfo.ymid             = center(2); % vertical center
ScreenInfo.numPixels_perCM  = 7.5;
ScreenInfo.liftingYaxis     = 300; 

%fixation locations
ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;

%% initialise serial object
%motorArduino = serial('/dev/cu.usbmodem14141'); 
motorArduino = serial('/dev/cu.usbmodemFA1341');
motorArduino.Baudrate=9600;
motorArduino.StopBits=1;
motorArduino.Terminator='LF';
motorArduino.Parity='none';
motorArduino.FlowControl='none';
% open for usage
fopen(motorArduino);
% pause for process delay
pause(2);
%set number of steps to be moved
noOfSteps = 3200;  %how many microsteps
%how many steps do we need to move one centimeter

%% open loudspeakers and create sound stimuli 
PsychDefaultSetup(2);
% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(3).DeviceIndex;

%% create auditory stimuli
%white noise for mask noise
AudInfo.fs              = 44100;
audioSamples            = linspace(1,AudInfo.fs,AudInfo.fs);
maskNoiseDuration       = 2; %the duration of the mask sound depends the total steps
MotorNoiseRepeated      = MakeMaskingSound(AudInfo.fs*maskNoiseDuration);
%set sound duration for both sound stimulus and mask noise
duration_mask           = length(audioSamples)*maskNoiseDuration;
timeline_mask           = linspace(1,duration_mask,duration_mask);
%generate white noise (one audio output channel)
carrierSound_mask       = randn(1, max(timeline_mask)); 
%create adaptation sound (only one loudspeaker will play the white noise)
AudInfo.intensity_MNR   = 100; %how loud do we want the recorded motor noise be
AudInfo.MaskNoise       = [carrierSound_mask+AudInfo.intensity_MNR.*MotorNoiseRepeated;
                            zeros(size(carrierSound_mask))]; 
%windowed: sineWindow_mask.*carrierSound_mask %non-windowed: carrierSound_mask

%for gaussion white noise
standardFrequency_gwn       = 10;
AudInfo.adaptationDuration  = 0.1; %0.05 %the burst of sound will be displayed for 40 milliseconds
duration_gwn                = length(audioSamples)*AudInfo.adaptationDuration;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs); 
carrierSound_gwn            = randn(1, max(timeline_gwn));
AudInfo.intensity_GWN       = 15;
AudInfo.GaussianWhiteNoise  = [zeros(size(carrierSound_gwn));...
                                 AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn]; 
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

%% Specify the experiment information
ExpInfo.sittingDistance = 105;
ExpInfo.numTrials       = 1; 
ExpInfo.testLocations   = 4;  %V = -10.5, -3.5, 3.5, 10.5
ExpInfo.numStaircases   = 8;  %2 for each test location

%% define visual stimuli
VSinfo.locations_deg    = -12:8:12;
VSinfo.locations_cm     = round(tan(deg2rad(VSinfo.locations_deg))*...
                            ExpInfo.sittingDistance,4);%[-13.8 -4.6 4.6 13.8]; %in cm
VSinfo.numFrames        = 6;   %the number of frames for the visual stimulus (100 ms) 
VSinfo.width            = 401; %(pixel) Increasing this value will make the cloud more blurry
VSinfo.boxSize          = 201; %This is the box size for each cloud.
VSinfo.intensity        = 10;  %This determines the height of the clouds. 
                               %Lowering this value will make them have lower contrast

%set the parameters for the visual stimuli, which consist of 10 gaussian blobs
VSinfo.blackScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
x                  = 1:1:VSinfo.boxSize; y = x;
[X,Y]              = meshgrid(x,y);
cloud              = 1e2.*mvnpdf([X(:) Y(:)],[median(x) median(y)],...
                        [VSinfo.width 0; 0 VSinfo.width]);
VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(cloud,length(x),length(y)); 
VSinfo.blk_texture = Screen('MakeTexture', windowPtr, VSinfo.blackScreen,[],[],[],2);  

%% auditory information
AudInfo.initialLoc_deg     = -7.5;
AudInfo.initialOffset_deg  = [1.25, 5.25, 12.75, 10.75; 10.75, 12.75, 5.25, 1.25];
AudInfo.StepSizes          = [1.5, 0.75, 0.375]; %in deg
AudInfo.waitTime           = 4;
%initialize the distance for all 4 test locations (8 staircases)
temp_distance       = repmat(VSinfo.locations_deg,[2 1]);
initial_offset      = AudInfo.initialOffset_deg(:).*repmat([-1 1],[1 ExpInfo.testLocations])';
AudInfo.Distance    = temp_distance(:) + initial_offset;

%% Practice trials
%record data for easy trials
ExpInfo.numEasyTrialsPerS = 1;
ExpInfo.numEasyTrials     = ExpInfo.numEasyTrialsPerS*ExpInfo.numStaircases;

%2nd row: the V loc in deg
%3rd row: the A loc in deg
%4th row: order (1: A first; 2: V first)
%5th row: response (-1: A is to the left of V; 1: A is to the right of V)
%6th row: Response time
D_easyTrials      = NaN(6,ExpInfo.numEasyTrials);
D_easyTrials(1,:) = randperm(ExpInfo.numStaircases, ExpInfo.numStaircases);
idx_testLoc       = fix((D_easyTrials(1,:)+1)/2);
D_easyTrials(2,:) = VSinfo.locations_deg(idx_testLoc);
D_easyTrials(3,:) = AudInfo.Distance(D_easyTrials(1,:),1);
D_easyTrials(4,:) = Shuffle(repmat([1 2],[1 ExpInfo.numEasyTrials/2]));

%% Run the experiment by calling the function InterleavedStaircase
%start the experiment
DrawFormattedText(windowPtr, 'Press any button to start the experiment.',...
	'center', ScreenInfo.yaxis - ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr); KbWait(-3); WaitSecs(1);
Screen('Flip',windowPtr);

currentAloc = AudInfo.initialLoc_deg; %in degs
for i = 1:ExpInfo.numEasyTrials
    [D_easyTrials(5,i),D_easyTrials(6,i)] = EasyTrials(currentAloc,...
        D_easyTrials(3,i),D_easyTrials(2,i),D_easyTrials(4,i),ScreenInfo,...
        VSinfo,AudInfo,ExpInfo,motorArduino,noOfSteps,pahandle,windowPtr); 
    currentAloc = D_easyTrials(3,i);
end

%% Move the arduino to the leftmost place
if abs(currentAloc - AudInfo.initialLoc_deg)>=0.01
    distance_goingBack = round((tan(deg2rad(AudInfo.initialLoc_deg))-...
        tan(deg2rad(currentAloc)))*ExpInfo.sittingDistance/3,2); 
    if distance_goingBack < 0 
        fprintf(motorArduino,['%c','%d'], ['p',noOfSteps*abs(distance_goingBack)]);
    else
        fprintf(motorArduino,['%c','%d'], ['n',noOfSteps*abs(distance_goingBack)]);
    end 
end
delete(motorArduino)

%% Save data and end the experiment
A_aligns_V_practice_data = {ExpInfo,ScreenInfo,VSinfo,AudInfo,D_easyTrials};
save(out1FileName,'A_aligns_V_practice_data');
Screen('CloseAll');