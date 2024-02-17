speakerSpan = 65.5 * 2; % cm
sittingDist = 113; % cm
screenWidth = 170; % cm
screenX = 1024; % pixel
screenMid = mean(1:screenX);
numSpeakerIntervals = 15; % 15 intervals between 16 speakers
cmPerAudInd = speakerSpan / numSpeakerIntervals;
pixelPerCm = screenX / screenWidth;

ExpInfo.audLevel = [5 7 10 12];%[5,7,10,12];
mostLRUsed = [min(ExpInfo.audLevel),max(ExpInfo.audLevel)];
angleRange = rad2deg(atan((mostLRUsed(2) - mostLRUsed(1)) / 2 * cmPerAudInd / sittingDist));
angleRange = [-angleRange, angleRange];

sA = round((ExpInfo.audLevel - 8.5) .* cmPerAudInd .* pixelPerCm + screenMid); % in pixel. 8.5 being middle point
sV = sA;
sAV = combvec(sA, sV);
num_s = size(sAV,2); % number of conditions
nT = 1000;

aA        = 1;
bA        = 0;

sigA  = 100;
sigV = 25;

muP = screenMid;
sigP = 1000; % arbitarily using screen width here

pCommon = 0.57; % only 1/4 of the trials are common cause so I assume this here

ds        = {'model averaging', 'model selection', 'probability matching'};
nDS       = length(ds);

bimResp = NaN(num_s, length(ds), 3, nT);
% num_s      = stimulus combination
% length(ds) = 3 decision strategies (1. averaging, 2. selection, 3. matching)
% 3          = 2 modalities (1. auditory, 2. visual) and post_common
% nT         = num_trial
x = 1:1024; % the screen pixel space
%%
% for debug
% sA = sAV(1,i);
% sV = sAV(2,i);

for i = 1:num_s
    [bimResp(i,:,:,:),pdf(i)] = simResp_CI_solutions_Multi(pCommon, nT, sAV(1,i),...
        sAV(2,i), aA, bA, sigA, sigV, muP, sigP,x);
end


%% Model Averaging could be using...
% separate cause pdf
% common cause pdf
% averaged posterior pdf
maxPoint = 1;
minPoint = 0.01;
elbow = 256; % 1/4 of ScreenX is what we are using right now
strategy = 1;
for i = 1:num_s
radiusMA(i).A_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusMA(i).A_C2 = eGain(pdf(i).sHat_A_C2, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusMA(i).A_AvgPost = eGain(pdf(i).Aaveraging, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);

radiusMA(i).V_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusMA(i).V_C2 = eGain(pdf(i).sHat_V_C2, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusMA(i).V_AvgPost = eGain(pdf(i).Vaveraging, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
end

%% Model Selection could be using...
% separate cause pdf
% common cause pdf
% the selected pdf
maxPoint = 1;
minPoint = 0.01;
elbow = 256; % 1/4 of ScreenX is what we are using right now
strategy = 2;
for i = 1:num_s
radiusMS(i).A_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusMS(i).A_C2 = eGain(pdf(i).sHat_A_C2, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusMS(i).A_MSselect = eGain(pdf(i).A_selected, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);

radiusMS(i).V_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusMS(i).V_C2 = eGain(pdf(i).sHat_V_C2, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusMS(i).V_MSselect = eGain(pdf(i).V_selected, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
end

%% Probability Matching could be using...
% separate cause pdf
% common cause pdf
% selected pdf
% the mixture pdf
maxPoint = 1;
minPoint = 0.01;
elbow = 256; % 1/4 of ScreenX is what we are using right now
strategy = 3;
for i = 1:num_s
radiusPM(i).A_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).A_C2 = eGain(pdf(i).sHat_A_C2, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).A_PMselect = eGain(pdf(i).A_PMselected, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).A_PMmixture = eGain(pdf(i).AMatching, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);

radiusPM(i).V_C1 = eGain(pdf(i).sHat_C1, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).V_C2 = eGain(pdf(i).sHat_V_C2, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).V_PMselect = eGain(pdf(i).V_PMselected, round(squeeze(bimResp(i,strategy,2,:))), maxPoint, minPoint, elbow, screenX);
radiusPM(i).V_PMmixture = eGain(pdf(i).VMatching, round(squeeze(bimResp(i,strategy,1,:))), maxPoint, minPoint, elbow, screenX);
end

%% Plot simulated data
figure(1)
cond_indice = reshape(1:num_s,4,4);
modality = 1;
for i = 1:4
    for j = 1:4
        cond_ind = cond_indice(j,i);
        subplot(4,4,cond_ind)
        % h1 = histogram(squeeze(bimResp(cond_ind,1,modality,:)),'FaceAlpha',0.6);
        hold on
        h2 = histogram(squeeze(bimResp(cond_ind,2,modality,:)),'FaceAlpha',0.6);
        h3 = histogram(squeeze(bimResp(cond_ind,3,modality,:)),'FaceAlpha',0.6);
        xlim([100,900])
    end
end
hold off
%% confidence radius
figure(2)
cond_indice = reshape(1:num_s,4,4);
plotInd = 1;
for i = 1:4
    for j = 1:4
        cond_ind = cond_indice(j,i);
        subplot(4,4,plotInd)
        % h1 = histogram(squeeze(radiusMA(cond_ind).A_AvgPost),'FaceAlpha',0.6);
        hold on
        h2 = histogram(squeeze(radiusMS(cond_ind).A_MSselect),'FaceAlpha',0.6);
        h3 = histogram(squeeze(radiusPM(cond_ind).A_PMselect),'FaceAlpha',0.6);
        % xlim([20,50])
        title(['sAV = ' num2str(sAV(:,cond_ind)')])

        plotInd = plotInd + 1;
    end
end
%% posterior of common cause
figure(3)
cond_indice = reshape(1:num_s,4,4);
modality = 2;
plotInd = 1;
for i = 1:4
    for j = 1:4
        cond_ind = cond_indice(j,i)
        % cond_ind = plotInd;
        subplot(4,4,plotInd)
        h1 = histogram(squeeze(bimResp(cond_ind,1,3,:)),'FaceAlpha',0.6);
        % xlim([20,50])
        title(['sAV = ' num2str(sAV(:,cond_ind)')])
        xlabel('Post_{common}')
        ylabel('Count')
        plotInd = plotInd + 1;
    end
end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
