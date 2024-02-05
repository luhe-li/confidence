speakerSpan = 65.5 * 2; % cm
sittingDist = 113; % cm
screenWidth = 170; % cm
screenX = 1024; % pixel
screenMid = mean(1:screenX);
numSpeakerIntervals = 15; % 15 intervals between 16 speakers
cmPerAudInd = speakerSpan / numSpeakerIntervals;
pixelPerCm = screenX / screenWidth;

ExpInfo.audLevel = [5,7,10,12];
mostLRUsed = [min(ExpInfo.audLevel),max(ExpInfo.audLevel)];
angleRange = rad2deg(atan((mostLRUsed(2) - mostLRUsed(1)) / 2 * cmPerAudInd / sittingDist));
angleRange = [-angleRange, angleRange];

sA = round((ExpInfo.audLevel - 8.5) .* cmPerAudInd .* pixelPerCm + screenMid); % in pixel 
sV = sA;
sAV = combvec(sA, sV);
num_s = size(sAV,2);
nT = 100;

aA        = 1;
bA        = 0;

sigA  = 80;
sigV = 20;

muP = screenMid;
sigP = 170; % arbitarily using screen width here

pCommon = 0.25; % only 1/4 of the trials are common cause so I assume this here

ds        = {'model averaging', 'model selection', 'probability matching'};
nDS       = length(ds);

bimResp = NaN(num_s, length(ds), 2, nT);
x = 1:1024;
%%
for i = 1:num_s

    % stimulus combination
    % 3 decision strategies (1. averaging, 2. selection, 3. matching)
    % 2 modalities (1. auditory, 2. visual)
    % num_trial
    bimResp(i,:,:,:) = simResp_CI_solutions_Multi(pCommon, nT, sAV(1,i),...
        sAV(2,i), aA, bA, sigA, sigV, muP, sigP);

end






