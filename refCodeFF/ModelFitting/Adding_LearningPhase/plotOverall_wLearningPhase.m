%% load data
clear all; close all; clc
global subjIs
subjNs     = [3,4,5,6,8,9,11,12,13,15,16,17,18,19,20];
subjIs     = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL','RE','MM'};
modelTypes = {'','_diffSigma'};
strategy   = {'strategyMAP_MA_strategyUnity_posteriorC1',...
              'strategyMAP_MA_strategyUnity_measurements',...
              'strategyMAP_MA_strategyUnity_MAP',...
              'strategyMAP_MS_strategyUnity_posteriorC1',...
              'strategyMAP_MS_strategyUnity_measurements'};
cond       = {'cong', 'incong'};
phase      = {'pre', 'post'};
modality   = {'Auditory','Visual'};
idx_mat    = {[8,9;10,11];[12,13;14,15]};

%define useful variables and initialize matrices
global lenC lenP lenM lenMT lenDS nS
lenC               = length(cond); 
lenP               = length(phase); 
lenM               = length(modality);
lenMT              = length(modelTypes);
lenDS              = length(strategy);
nS                 = length(subjNs);
[minNLL, deltaNLL] = deal(NaN(nS, lenMT, lenDS));
estimatedP         = cell(1,nS);
pC1_all            = NaN(nS, lenMT, lenDS, lenC, lenP);
for i = 1:nS   
    %add folder path
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'ModelFitting/Fits/', subjIs{i},'/Fits_wLearningPhase']));
    for j = 1:lenMT
        for k = 1:lenDS
            %load data
            C = load(['ModelFitting_updatePrior_fullModel_', strategy{k},...
                '_wLearningPhase', modelTypes{j}, '_sub', num2str(subjNs(i)),...
                '.mat'], 'ModelFitting');
            estimatedP{i}{j,k} = C.ModelFitting{end}.P_f;
            minNLL(i,j,k)      = C.ModelFitting{end}.minNLL_f;
            pC1_all(i,j,k,1,:) = estimatedP{i}{j,k}(idx_mat{j}(1,:)); %cong condition
            pC1_all(i,j,k,2,:) = estimatedP{i}{j,k}(idx_mat{j}(2,:)); %incong condition
        end
    end
    deltaNLL(i,:,:) = minNLL(i,:,:) - min(min(minNLL(i,:,:)));
end

%% plot AIC
green           = [76,153,0]./255;
colorMap        = [linspace(green(1),1,255)', linspace(green(2),1,255)',...
                    linspace(green(3),1,255)'];
getAIC          = @(nLL, k) 2.*nLL+ 2.*k;
numParams       = [13,14,14,13,14; 17,18,18,17,18]; 
%numTotalTrials  = 320 + 240 + 320*2*2; 
%initialize
[AIC, deltaAIC]         = deal(NaN(nS, lenMT, lenDS));
subjI_bestM             = cell(1,nS); %best model name
subjI_bestM_sub         = NaN(nS, 2); %subscripts
subjI_bestS_saturated   = NaN(nS);
%the pCommon, sigma_AV_A, sigma_AV_V given the best-fitting decision strategy
pC1_saturated           = NaN(nS, lenC, lenP); 
sigma_saturated         = NaN(nS, lenM, lenP + 1);
sigma_diff_saturated    = NaN(nS, lenM, lenC);
bool_plotAIC            = 0;
for i = 1:nS
    %calculate delta AIC for each subject
    AIC(i,:,:)           = getAIC(squeeze(minNLL(i,:,:)), numParams);
    deltaAIC(i,:,:)      = AIC(i,:,:) - min(min(AIC(i,:,:)));
    %find the best model (the model that corresponds to 0 delta AIC)
    idx_bestM            = find(squeeze(deltaAIC(i,:,:)) == 0);
    [row, col]           = ind2sub([lenMT, lenDS], idx_bestM);
    subjI_bestM{i}       = [subjIs{i},modelTypes{row},'-', strategy{col}];
    subjI_bestM_sub(i,:) = [row, col];
    disp(subjI_bestM{i}); %display it
    
    %given saturated model, find the best strategy
    [~,idx_bestS]               = min(deltaAIC(i,2,:));
    subjI_bestS_saturated(i)    = idx_bestS;
    %get the sigma for all participants
    pC1_saturated(i,:,:)        = squeeze(pC1_all(i,2, idx_bestS,:,:));
    sigma_saturated(i,1,:)      = estimatedP{i}{2, idx_bestS}(5:7); %A
    sigma_saturated(i,2,:)      = estimatedP{i}{2, idx_bestS}(8:10); %V
    sigma_diff_saturated(i,:,1) = squeeze((sigma_saturated(i,:,2) - sigma_saturated(i,:,1))); %cong
    sigma_diff_saturated(i,:,2) = squeeze((sigma_saturated(i,:,end) - sigma_saturated(i,:,1))); %incong
    %plot
    if bool_plotAIC == 1
        figure
        heatmap(round(squeeze(deltaAIC(i,:,:)),1),'ColorbarVisible', 'on',...
            'YLabel', sprintf('Assumption \nabout sigmas'),...
            'XLabel', sprintf('Assumption about strategies'), 'YData',...
            {'no change', 'changes'}, 'XData', ...
            {sprintf('MA-posterior_{C=1}'), sprintf('MA-measurements'),...
            sprintf('MA-estimates'),sprintf('MS-posterior_{C=1}'),...
            sprintf('MS-measurements')},'GridVisible','off','Colormap',colorMap); 
        caxis([0, 10]); set(gca,'FontSize',18); 
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.57, 0.3]);
        set(gcf,'PaperUnits','centimeters','PaperSize',[40 10]);
        %saveas(gcf, ['ModelComparison_deltaAIC_',subjIs{i}], 'pdf'); 
    end
end
pC1_diff_saturated = squeeze(pC1_saturated(:,:,end) - pC1_saturated(:,:,1));
%size: nS x lenMT x 2 lenC

%% count best-fit model
AIC_thres   = 1e-4;
counts_M    = squeeze(sum(deltaAIC < AIC_thres,1));
org         = [255,140,0]./255;
colorMap    = [linspace(1,org(1),255)',linspace(1,org(2),255)',linspace(1,org(3),255)'];
lbls        = {sprintf('MA-posterior_{C=1}'), sprintf('MA-measurements'),...
               sprintf('MA-estimates'),sprintf('MS-posterior_{C=1}'),...
               sprintf('MS-measurements')};
           
figure
heatmap(counts_M,'ColorbarVisible', 'on','YLabel', ...
    sprintf('Assumption \nabout sigmas'),'XLabel', ...
    sprintf('Assumption about strategies'), 'YData',...
    {'no change', 'changes'}, 'XData', lbls,'GridVisible','off',...
    'Colormap',colorMap); 
caxis([0, nS]); set(gca,'FontSize',18); %25
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.57, 0.3]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 10]);
% saveas(gcf, ['ModelComparison_counts_bestModel_DeltaAICthres',num2str(AIC_thres)], 'pdf'); 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.60, 0.32]);
% set(gcf,'PaperUnits','centimeters','PaperSize',[45 13]);
% saveas(gcf, ['ModelComparison_counts_bestModel_DeltaAICthres',...
%     num2str(AIC_thres),'_poster'], 'pdf'); 

%% plotting the change of sigma and pCommon
plot_bestfitPC1_diffSigma(pC1_diff_saturated, sigma_diff_saturated,...
    deltaAIC, subjI_bestS_saturated, 1)

