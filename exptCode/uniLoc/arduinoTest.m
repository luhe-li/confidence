%% Enter experiment info
clear; close all;  rng('Shuffle');

Arduino = serial('/dev/cu.usbmodem14101','BaudRate',115200);
fopen(Arduino);

%% Auditory set up

% open speakers and create sound stimuli
PsychDefaultSetup(2);

% get correct sound card
InitializePsychSound
devices = PsychPortAudio('GetDevices');
our_device=devices(end).DeviceIndex;

%%
% Gaussian white noise
AudInfo.fs                  = 44100;
audioSamples                = linspace(1,AudInfo.fs,AudInfo.fs);
standardFrequency_gwn       = 100;
AudInfo.stimDura            = 5; % sec
duration_gwn                = length(audioSamples)*AudInfo.stimDura;
timeline_gwn                = linspace(1,duration_gwn,duration_gwn);
sineWindow_gwn              = sin(standardFrequency_gwn/2*2*pi*timeline_gwn/AudInfo.fs);
carrierSound_gwn            = randn(1, numel(timeline_gwn));
AudInfo.intensity_GWN       = 1; % too loud for debugging, originally 15
% AudInfo.GaussianWhiteNoise  = [zeros(size(carrierSound_gwn)); zeros(size(carrierSound_gwn))];
AudInfo.GaussianWhiteNoise  = [AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn; AudInfo.intensity_GWN.*sineWindow_gwn.*carrierSound_gwn];
pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device

aa = repmat([6,11],[1,20]);

%% audio test
for i = aa
    testSpeaker = i;
    input_on = ['<',num2str(1),':',num2str(i),'>']; %arduino takes input in this format
    fprintf(Arduino,input_on);
    PsychPortAudio('FillBuffer',pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('Start',pahandle,1,0,0);
    WaitSecs(5.5);
    input_off = ['<',num2str(0),':',num2str(i),'>'];
    fprintf(Arduino,input_off);
    PsychPortAudio('Stop',pahandle);
    WaitSecs(1);
end