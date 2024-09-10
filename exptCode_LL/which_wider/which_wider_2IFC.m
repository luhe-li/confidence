% "Which one is wider?" 2IFC task

%% experiment info

clear; close all; clc; rng('Shuffle');

ExpInfo.sub_init = [];
while isempty(ExpInfo.sub_init) == 1
    try ExpInfo.sub_init = input('Participant Initial: ','s') ;
        ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
    catch
    end
end

switch ExpInfo.practice
    case 1
        out_finm = sprintf('wide_sub%s_ses-%s', ExpInfo.sub_init);
    case 2
        out_finm = sprintf('wide_practice_sub%s_ses-%s', ExpInfo.sub_init);
end

% path control
[project_dir, ~]  = fileparts(fileparts(pwd));
[git_dir, ~] = fileparts(project_dir);
out_dir = fullfile(project_dir, 'data', 'which_wider');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(git_dir, 'Psychtoolbox-3')))

% avoid rewriting data
if exist(fullfile(out_dir, [out_finm '.mat']), 'file')
    resp = input('Replace the existing file? Y/N', 's');
    if ~strcmp(resp,'Y')
        disp('Experiment terminated.')
        return
    end
end

% switch between debug mode
ExpInfo.mode  = 2; %input('Experiment mode: 1; Debug mode: 2#: ');
switch ExpInfo.mode
    case 1 % experiment mode
        windowSize = [];
        opacity = 1;
        HideCursor();
    case 2 % debug mode
        windowSize = [100 100 1000 600]; % open a smaller window
        opacity = 0.4;
end

%% screen Setup

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
[windowPtr,rect] = Screen('OpenWindow', screenNumber, [0,0,0], windowSize);
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

%% experiment setup

% expt info
ExpInfo.sig_dot = 4; % cm, sig_sample tested: 4, 8, 12, 16
ExpInfo.sig_sample = [0, 16];
ExpInfo.n_trial = 36; % for each staircase
ExpInfo.n_staircase = 2; % one up and one down
ExpInfo.n_block = 4;
ExpInfo.n_total_trial = ExpInfo.n_staircase * ExpInfo.n_trial;

% width/distance info
ExpInfo.speaker_cm = 65.5; % cm, left to center
ExpInfo.sitting_dist = 113; % cm
ExpInfo.screen_cm = 170; % cm
ExpInfo.px_per_cm = rect(3)/ExpInfo.screen_cm;

% choose speaker/visual location
ExpInfo.speaker_idx = 6:11;
ExpInfo.speaker_level_cm = linspace(-ExpInfo.speaker_cm, ExpInfo.speaker_cm, 16);
ExpInfo.sA_cm = ExpInfo.speaker_level_cm(ExpInfo.speaker_idx);
ExpInfo.sV_cm = ExpInfo.sA;

% duration
ExpInfo.n_frame = 12;
ExpInfo.t_fixation = 0.5;
ExpInfo.t_blank1 = 0.2;
ExpInfo.ITI = 0.3;
ExpInfo.tIFI = ScreenInfo.ifi;

%% define visual stimuli

VSinfo.width            = 401; %(pixel) Increasing this value will make the cloud more blurry
VSinfo.boxSize          = 201; %This is the box size for each cloud.
VSinfo.intensity        = 10;  %This determines the height of the clouds.
%Lowering this value will make them have lower contrast

%set the parameters for the visual stimuli
VSinfo.pblack              = 1/8; % set contrast to 1*1/8 for the "black" background, so it's not too dark and the projector doesn't complain
VSinfo.greyScreen          = VSinfo.pblack * ones(ScreenInfo.xaxis,ScreenInfo.yaxis)*255;
VSinfo.grey_texture        = Screen('MakeTexture', windowPtr, VSinfo.greyScreen,[],[],[],2);
VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
x                  = 1:1:VSinfo.boxSize; y = x;
[X,Y]              = meshgrid(x,y);
cloud              = 1e2.*mvnpdf([X(:) Y(:)],[median(x) median(y)],...
    [VSinfo.width 0; 0 VSinfo.width]);
% VSinfo.Cloud       = (1-VSinfo.pblack).*255.*VSinfo.intensity.*reshape(cloud,length(x),length(y));

%% initializing vectors that store the data

[Conditions,Response,Order,RTs] = deal(NaN(ExpInfo.n_staircase, ExpInfo.n_trial));
%------------------------------Conditions----------------------------------
%odd number: \sig_sample starts of 0, smaller than \sig_dot(1-up-2-down)
%even number: \sig_sample starts from the max tested value, larger than \sig_dot (2-up-1-down)
for i = 1:ExpInfo.n_trial
    Conditions(:,i) = randperm(ExpInfo.n_staircase, ExpInfo.n_staircase);
end
%---------------------------------Order------------------------------------
%1: present standard stimulus first
%2: present comparison stimulus first
for i = 1:ExpInfo.n_trial
    Order(:,i) = Shuffle(repmat([1,2],[1, ExpInfo.n_staircase/2]));
end
%-------------------------------Response-----------------------------------
%1: participants think the first is wider
%2: participants think the second stimulus is wider

%% Add easy trials to calculate lapse rate

% catch trials are evenly spread across blocks
ExpInfo.n_easy_trial_per_s  = 4;
trial_slc  = [];
for i = 1:ExpInfo.n_staircase
    trialIdx_staircase = reshape(i:ExpInfo.n_staircase:ExpInfo.n_total_trial, ...
        [ExpInfo.n_trial/ExpInfo.n_block, ExpInfo.n_block]);
    for j = 1:ExpInfo.n_easy_trial_per_s
        trial_slc = [trial_slc, trialIdx_staircase(randi(...
            ExpInfo.n_trial/ExpInfo.n_block),j)];
    end
end
%To see whether the catch trials are evenly spread out. Do:
T = zeros(ExpInfo.n_staircase, ExpInfo.n_trial);
T(trial_slc) = 1; imagesc(T);

%% run staircase

instruction = ['In the following session, you will see two visual stimuli in sequence on each trial.','\nPress 1 if you think the first stimulus is wider.', '\nPress 2 if you think the second stimulus is wider.'];

c                   = clock;
ExpInfo.start       = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
DrawFormattedText(windowPtr, instruction ,...
    'center',ScreenInfo.yaxis-500,[255 255 255]);
Screen('Flip',windowPtr);
KbWait(-3);
WaitSecs(1);

for i = 1:ExpInfo.n_trial

    for j = 1:ExpInfo.n_staircase

        idx_s = Conditions(j,i);
        if i == 1; history = []; else history = Response(idx_s,1:(i-1)); end

        [Response(idx_s,i),RTs(idx_s,i),Distance(idx_s,i+1)] = ...
            interleaved_staircase(idx_s,i,currentAloc,history,...
            Distance(idx_s,i),ScreenInfo,ExpInfo,VSinfo,AudInfo,...
            Order(idx_s,i),motorArduino,noOfSteps,pahandle,windowPtr);

        %inserted an easy trial for calculating the lapse rate later
        idx_t = (i-1)*ExpInfo.n_staircase+j;

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