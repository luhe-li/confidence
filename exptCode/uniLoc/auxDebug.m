%% Enter experiment info
clear; close all;  rng('Shuffle');

% Arduino = serial('/dev/cu.usbmodem14101','BaudRate',115200);
% fopen(Arduino);

%% Auditory set up

% open speakers and create sound stimuli
PsychDefaultSetup(2);

% get correct sound card
InitializePsychSound;
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn       = 20;
AudInfo.stimDura            = 10; % sec
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN       = 1; % too loud for debugging, originally 15
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn;
    zeros(1,duration_gwn)];
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

% aa = reshape([6:11; 6:11;6:11],1,[]);
aa = [repmat([1 2 15 16],[1, 40])];
%% audio test
for i = aa
        PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(10);
        PsychPortAudio('Stop',pahandle);
    WaitSecs(2 );
end