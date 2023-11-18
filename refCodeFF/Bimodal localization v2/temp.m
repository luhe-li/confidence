%% load data
clear all; close all; clc
addpath(genpath('/e/3.3/p3/hong/Desktop/Project2/Bimodal localization/data'));
subjN = 11;
cond  = {'congruent','incongruent'};
phase = {'pre','post'};
numTotalTrials = 320;
[AVpairs_order, unityResp, locResp] = deal(NaN(length(cond), length(phase),...
    numTotalTrials));
for i = 1:length(cond)
    for j = 1:length(phase)
        C = load(['BimodalLocalization_', phase{j},'_sub', num2str(subjN),...
            '_session',num2str(i),'.mat'],['BimodalLocalization_',phase{j},...
            '_data']);
        AVpairs_order(i,j,:) = eval(['C.BimodalLocalization_',phase{j},...
            '_data{1}.AVpairs_order']);
        unityResp(i,j,:) = eval(['C.BimodalLocalization_',phase{j},...
            '_data{end}.unity']);
        locResp(i,j,:) = eval(['C.BimodalLocalization_',phase{j},...
            '_data{end}.localization']);
    end
end
%other useful info
numTrialsPerPair    = C.BimodalLocalization_post_data{1}.numTrials;
A_loc              = C.BimodalLocalization_post_data{4}.Distance;
V_loc              = C.BimodalLocalization_post_data{3}.Distance;
AVpairs_allComb    = C.BimodalLocalization_post_data{1}.AVpairs_allComb;
                        %1st row: A locations; 2nd row: V locations

%% organize data
[unityMat,locRespMat] = deal(NaN(length(cond), length(phase), length(A_loc),...
                        length(V_loc),numTrialsPerPair));
[pC1_resp,locRespMean,locRespSTD] = deal(NaN(length(cond), length(phase),...
                        length(A_loc), length(V_loc)));
for i = 1:length(cond)
    for j = 1:length(phase)
        for k = 1:size(AVpairs_allComb,2)
            idx_samePair = (AVpairs_order(i,j,:) == k);
            unityResp_samePair = unityResp(i,j,idx_samePair);
            unityMat(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k),:) = ...
                unityResp_samePair;
            pC1_resp(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k)) = ...
                sum(unityResp_samePair==1)/numTrialsPerPair;
        end
    end
end

%% show figures
figure(1)
for i = 1:length(cond)
    for k = 1:length(V_loc)
        subplot(length(cond), length(V_loc), (i-1)*length(V_loc)+k);
        plot(A_loc, squeeze(pC1_resp(i,1,:,k)),'--','lineWidth',2, 'Color',...
            [135,206,235]./255); hold on
        plot(A_loc, squeeze(pC1_resp(i,2,:,k)),'-','lineWidth',2,'Color',...
            [65,105,225]./255); hold on
        box off; xticks(A_loc); xlim([A_loc(1), A_loc(end)]); xticks(A_loc);
        if k == 1; yticks(0:0.5:1); ylabel('P(reporting C=1)');
        else; yticks([]);end
        if i == 2; xlabel('A location (deg)'); end
        if i == 1 && k == 1; legend(phase); legend boxoff;end
        title(['V = ', num2str(V_loc(k)),char(176)]);
        set(gca,'FontSize',15);
    end
end








