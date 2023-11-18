%% Test all speakers to make sure everything works
clear all; close all;
addpath(genpath('/e/3.3/p3/hong/Desktop/Project5/Psychtoolbox'));
%% Initialise serial object
Arduino = serial('/dev/cu.usbmodemFD131','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
%open for usage
fopen(Arduino); 

%% Open speakers and create sound stimuli 
PsychDefaultSetup(2);
% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(3).DeviceIndex;

%% make auditory stimuli (beep)
% AudInfo.fs                           = 44100;
% AudInfo.stimDura                     = 0.1; 
% AudInfo.tf                           = 400; 
% AudInfo.intensity                    = 0.65;
% beep                                 = MakeBeep(AudInfo.tf, AudInfo.stimDura, AudInfo.fs);
% AudInfo.inBetweenBeep                = [AudInfo.intensity*beep; AudInfo.intensity*beep];
% AudInfo.Beep                         = [beep; beep];
% pahandle = PsychPortAudio('Open', our_device, [], [], [], 2); %open device

%% Gaussian white noise
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

%% Make visual stimuli (gaussian blob)
screens = Screen('Screens');
screenNumber = max(screens);
black = BlackIndex(screenNumber);
opacity = 1;
PsychDebugWindowConfiguration([], opacity)

[windowPtr, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
[center(1), center(2)]     = RectCenter(windowRect);
ScreenInfo.xmid            = center(1); % horizontal center
ScreenInfo.ymid            = center(2); % vertical center
ScreenInfo.numPixels_perCM = 6.2;
ScreenInfo.liftingYaxis    = 270; 
ifi = Screen('GetFlipInterval', windowPtr);
waitframes = 1;

%fixation locations
ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;

%calculate visual angle
ExpInfo.sittingDistance              = 113.0;
ExpInfo.leftspeaker2center           = 65.5;
ExpInfo.rightspeaker2center          = 65.5;
ExpInfo.leftmostVisualAngle          = (180/pi) * atan(ExpInfo.leftspeaker2center / ...
                                       ExpInfo.sittingDistance);
ExpInfo.rightmostVisualAngle         = (180/pi) * atan(ExpInfo.leftspeaker2center / ...
                                       ExpInfo.sittingDistance);
                                   
% intensity of standard stimulus
VSinfo.scaling                       = 0.4; % a ratio between 0 to 1 to be multipled by 255
pblack                               = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
% define the visual stimuli
VSinfo.Distance                      = linspace(-30,31,16); %in deg
VSinfo.numLocs                       = length(VSinfo.Distance);
VSinfo.numFrames                     = 6;
VSinfo.duration                      = VSinfo.numFrames * ifi;%s
VSinfo.width                         = 201; %(pixel) Increasing this value will make the cloud more blurry (arbituary value)
VSinfo.boxSize                       = 101; %This is the box size for each cloud (arbituary value)
%set the parameters for the visual stimuli
VSinfo.blackBackground               = pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis);
VSinfo.transCanvas                   = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
x                                    = 1:1:VSinfo.boxSize; y = x;
VSinfo.x                             = x; VSinfo.y = y;
[X,Y]                                = meshgrid(x,y);
cloud_temp                           = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
                                       [VSinfo.width 0; 0 VSinfo.width]);
pscale                               = (1-pblack)/max(cloud_temp); % the max contrast of the blob adds the background contrast should <= 1
cloud_temp                           = cloud_temp .* pscale;
VSinfo.Cloud                         = VSinfo.scaling.*reshape(cloud_temp,length(x),length(y));
VSinfo.blk_texture                   = Screen('MakeTexture', windowPtr, VSinfo.blackBackground,[],[],[],2);
                                   
                                   
%% test all speakers (present AV pairs)
nTrial = 2; %5;
nEvents = 31; %5;

% test stimuli 1-16 in order
left2right = 1:1:31;
right2left = 31:-1:1;
audTrain = [left2right; right2left];
visTrain = [left2right; right2left];

% % test stimuli 1-5 rand order no discrepancy 
% allSeqs = perms(1:5);
% audTrain = allSeqs(1:5,:);
% visTrain = allSeqs(6:10,:);

% % test rand order with discrepancy
% allSeqs = perms(1:5);
% audTrain = allSeqs(1:5,:);
% visTrain = allSeqs(6:10,:)+7; % +7; +11


for i=1:nTrial
    % fixation
    Screen('FillRect', windowPtr,[255 0 0], [ScreenInfo.x1_lb,...
        ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
    Screen('FillRect', windowPtr,[255 0 0], [ScreenInfo.x2_lb,...
        ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
    Screen('Flip',windowPtr); WaitSecs(1);

    % blank screen
    Screen('Flip',windowPtr);
    Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    WaitSecs(0.5);
    
    for j=1:nEvents
    %Calculate the coordinates of the target stimuli
        VSinfo.arrangedLocs_deg = VSinfo.Distance(visTrain(i,j));
        VSinfo.arrangedLocs_cm  = round(tan(deg2rad(VSinfo.arrangedLocs_deg)).*ExpInfo.sittingDistance,2);
        targetLoc = round(ScreenInfo.xmid + ScreenInfo.numPixels_perCM.*...
                    VSinfo.arrangedLocs_cm);
        %Make visual stimuli
        blob_coordinates = [targetLoc, ScreenInfo.liftingYaxis];    
        dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

        %----------------------------------------------------------------------
        %---------------------display audiovisual stimuli----------------------
        %----------------------------------------------------------------------

        if rem(audTrain(i,j),2) == 1 %turn on one speaker, full volume
           %present the audiovisual single event pair   
            input_on = ['<',num2str(1),':',num2str((audTrain(i,j)+1)/2),'>']; %arduino takes input in this format
            fprintf(Arduino,input_on);
            PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
            PsychPortAudio('Start',pahandle,1,0,0);

                for kk = 1:VSinfo.numFrames 
                        Screen('DrawTexture',windowPtr,dotCloud,[],...
                            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
                        Screen('Flip',windowPtr);
                end 
                for k = 1:12 %300ms blank screen
                        Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
                        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
                        Screen('Flip',windowPtr);
                end

            WaitSecs(0.10)
            input_off = ['<',num2str(0),':',num2str((audTrain(i,j)+1)/2),'>'];
            fprintf(Arduino,input_off);   
            PsychPortAudio('Stop',pahandle);
            WaitSecs(0.10)

        else %play sound in-between two speakers, half volume
            %present the audiovisual single event pair    
            input_on = ['<',num2str(1),':',num2str(audTrain(i,j)/2),',',...
                num2str(audTrain(i,j)/2+1),'>'];
            fprintf(Arduino,input_on);
            PsychPortAudio('FillBuffer',pahandle, AudInfo.inBetweenGWN);
            PsychPortAudio('Start',pahandle,1,0,0);

                for kk = 1:VSinfo.numFrames 
                        Screen('DrawTexture',windowPtr,dotCloud,[],...
                            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
                        Screen('Flip',windowPtr);
                end 

                for k = 1:12 %300ms blank screen
                        Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
                        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
                        Screen('Flip',windowPtr);
                end
            WaitSecs(0.10)
            input_off = ['<',num2str(0),':',num2str(audTrain(i,j)/2),',',...
                num2str(audTrain(i,j)/2+1),'>'];
            fprintf(Arduino,input_off);   
            PsychPortAudio('Stop',pahandle);
            WaitSecs(0.10)
        end
% %           present the audiovisual single event pair   
%             input_on = ['<',num2str(1),':',num2str((audTrain(i,j)+1)/2),'>'];
%             fprintf(Arduino,input_on);
%             PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
%             Screen('DrawTexture',windowPtr,dotCloud,[],...
%                             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%             Screen('Flip',windowPtr);
%             PsychPortAudio('Start',pahandle,1,0,0);  
% 
%                 for kk = 1:VSinfo.numFrames 
%                         Screen('DrawTexture',windowPtr,dotCloud,[],...
%                             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
%                         Screen('Flip',windowPtr);
%                 end 
%                 for k = 1:6 %300ms blank screen
%                         Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
%                         [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
%                         Screen('Flip',windowPtr);
%                 end
% 
%                 WaitSecs(0.1);
%             Screen('DrawTexture',windowPtr, VSinfo.blk_texture,[],...
%                         [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
%             Screen('Flip',windowPtr);
%             input_off = ['<',num2str(0),':',num2str((audTrain(i,j)+1)/2),'>'];
%             fprintf(Arduino,input_off);   
%             PsychPortAudio('Stop',pahandle);
%             WaitSecs(0.1)
          
    end 
   WaitSecs(2)
end

fclose(Arduino);
delete(Arduino)
Screen('CloseAll');
