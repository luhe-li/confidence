%% Enter experiment info
clear; close all;  rng('Shuffle');

ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Participant ID#: ') ;
        ExpInfo.session = input('Session: A/V1/V2#: ','s');
        ExpInfo.practice  = input('Main expt: 1; Practice: 2#: ');
    catch
    end
end

switch ExpInfo.practice
    case 1
        outFileName = sprintf('uniLoc_sub%i_ses-%s', ExpInfo.subjID, ExpInfo.session);
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
        if strcmp(ExpInfo.session, 'A')
        %Arduino = serial('/dev/cu.usbmodemFD131','BaudRate',115200); % make sure this value matches with the baudrate in the arduino code
        Arduino = serial('/dev/cu.usbmodem14101','BaudRate',115200);
        fopen(Arduino);
        end
    case 2 % debug mode
        windowSize = [100 100 1000 600]; % open a smaller window
        opacity = 0.4;
end
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
AudInfo.stimDura            = 0.033; % sec
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN       = 1; % too loud for debugging, originally 15
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;...
    AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

%% audio test
for i = 1:16
    testSpeaker = i;
    input_on = ['<',num2str(1),':',num2str(testSpeaker),'>']; %arduino takes input in this format
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(0.1);
    input_off = ['<',num2str(0),':',num2str(testSpeaker),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(0.1);
end