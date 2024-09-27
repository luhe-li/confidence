
%% Enter experiment info
clear; close all;  rng('Shuffle');

ExpInfo.subjInit = [];
while isempty(ExpInfo.subjInit) == 1
    try ExpInfo.subjInit = input('Participant Initial#: ','s') ;
        ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
    catch
    end
end

 switch ExpInfo.practice
     case 1
         outFileName = sprintf('point_sub-%s', ExpInfo.subjInit);
         ExpInfo.nRep = 20; % number of trial per condition level
         ExpInfo.numBlocks = 4;
     case 2
         outFileName = sprintf('point_practice_sub-%s', ExpInfo.subjInit);
         ExpInfo.nRep = 4; % number of trial per condition level
         ExpInfo.numBlocks = 2;
 end

 % path control
 curDir = pwd;
[projectDir, ~]  = fileparts(fileparts(curDir));
[git_dir, ~] = fileparts(projectDir);
addpath(genpath(fullfile(projectDir, 'exptCode','exptUtility')));
addpath(genpath(fullfile(git_dir, 'Psychtoolbox-3')))
outDir = fullfile(projectDir, 'data','pointingTask');
if ~exist(outDir,'dir') mkdir(outDir); end

% avoid rewriting data
if exist(fullfile(outDir, [outFileName '.mat']), 'file')
    resp = input('Replace the existing file? Y/N', 's');
    if ~strcmp(resp,'Y')
        disp('Experiment terminated.')
        return
    end
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
end
ExpInfo.randA = reshape(ExpInfo.randA, [], 1)'; 

% randomized auditory stimulus location in speaker index, use this as the
% reference of visual stimulus location in cm/pixel/VA
ExpInfo.randAudIdx = ExpInfo.audLevel(ExpInfo.randA);
ExpInfo.randVisIdx = ExpInfo.randAudIdx;

% location of speakers in CM, visual angle, and pixel
ExpInfo.sittingDistance              = 113.0; %cm
ExpInfo.numSpeaker                   = 16;
ExpInfo.LRmostSpeakers2center        = 65.5; %cm
ExpInfo.LRmostVisualAngle            = (180/pi) * atan(ExpInfo.LRmostSpeakers2center / ...
    ExpInfo.sittingDistance);
ExpInfo.speakerLocCM = linspace(-ExpInfo.LRmostSpeakers2center, ExpInfo.LRmostSpeakers2center, ExpInfo.numSpeaker);

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
ExpInfo.tFixation = 0.5;
ExpInfo.tBlank1 = 0.2;
ExpInfo.tStimFrame = 2;
ExpInfo.tStim = ExpInfo.tStimFrame * ScreenInfo.ifi;
ExpInfo.tIFI = ScreenInfo.ifi;
ExpInfo.tITI = 0.5;

% dial setup
ExpInfo.dialScaler = 2;
ExpInfo.conf_bar_height = 100;

%% make visual stimuli

% create the blob visual stimulus
VSinfo.stimFrame                     = ExpInfo.tStimFrame;
VSinfo.width                         = 31; %(pixel) Increasing this value will make the cloud more blurry
VSinfo.boxSize                       = 31; %This is the box size for each cloud.
x                                    = 1:1:VSinfo.boxSize; y = x;
VSinfo.x                             = x; VSinfo.y = y;
[X,Y]                                = meshgrid(x,y);
VSinfo.cloud                         = mvnpdf([X(:) Y(:)],[median(x) median(y)],[VSinfo.width 0; 0 VSinfo.width]);
VSinfo.pblack                        = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
pscale                               = (1-VSinfo.pblack)/max(VSinfo.cloud); % the max contrast of the blob adds the background contrast should <= 1
temp_cloud                           = VSinfo.cloud .* pscale;
VSinfo.Cloud                         = 255.*reshape(temp_cloud,length(x),length(y));
VSinfo.jitter_lb                     = -13*ScreenInfo.numPixels_perCM; % for jittering the visual target location, in pixel
VSinfo.jitter_ub                     = 13*ScreenInfo.numPixels_perCM;

% create background
VSinfo.greyScreen          = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture        = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen         = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
VSinfo.blackScreen = VSinfo.greyScreen;

%% Run the experiment

instruction = ['In the following session, you will be presented a white blob on each trial.','\nAfter the presentation, please use the cursor \nto locate the center of the blob.','\nUse the dial to adjust the net length, and press down the dial to confirm.','\nPress any key to start the pointing task.'];

%start the experiment
c                   = clock;
ExpInfo.start       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5) ,ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, instruction ,...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
WaitSecs(1);

for i = 1:ExpInfo.nTrials
    
    %% present stimuli
    
        SetMouse(ScreenInfo.xaxis*2, ScreenInfo.yaxis*2, windowPtr);
        HideCursor;
          Resp(i)= LocalizePointStim(i, ExpInfo,...
            ScreenInfo,VSinfo,windowPtr);    
        
    %% save by trial
    save(fullfile(outDir,outFileName),'Resp','ExpInfo','ScreenInfo','VSinfo');
    
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

% sort trials by location level
[~, temp_idx] = sort([Resp(1:end).target_idx]);
sortedResp = Resp(temp_idx); 
save(fullfile(outDir,outFileName),'Resp','sortedResp','ExpInfo','ScreenInfo','VSinfo');

%% display leaderoard

if  ExpInfo.practice == 1
    leaderboardText = updateLeaderboardPointing(outDir, ExpInfo, Resp);
else
    leaderboardText = 'This is the end of this session. Thank you!';
end
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, leaderboardText,...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
ShowCursor;
Screen('CloseAll');