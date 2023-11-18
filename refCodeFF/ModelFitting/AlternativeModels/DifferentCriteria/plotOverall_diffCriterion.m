%% load data
clear all; close all; clc
global subjIs
subjNs     = [3,4,5,6,8,9,11,12,13,15,16,17,18,19,20];
subjIs     = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL','RE','MM'};
modelTypes = {{'diffCriterion','4pC1_'},{'diffCriterion','2pC1_'},...
                {'fullModel',''},{'samePC1pre',''},{'samePC1cond',''}};
strategy   = {'strategyMAP_MA_strategyUnity_posteriorC1',...
              'strategyMAP_MA_strategyUnity_measurements',...
              'strategyMAP_MA_strategyUnity_MAP',...
              'strategyMAP_MS_strategyUnity_posteriorC1',...
              'strategyMAP_MS_strategyUnity_measurements'};
cond       = {'cong', 'incong'};
phase      = {'pre', 'post'};
modality   = {'Auditory','Visual'};
idx_pC1    = {[11,12;13,14],[11,11;12,12],[8,9;10,11],[8,9;8,10],[8,8;9,9]};
idx_crt    = {[7,8;9,10],[7,8;9,10],[7,7;7,7],[7,7;7,7],[7,7;7,7]};

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
[pC1_all, crt_all] = deal(NaN(nS, lenMT, lenDS, lenC, lenP));
for i = 1:nS   
    %add folder path
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'ModelFitting/Fits/', subjIs{i},'/Fits_alternativeModels/DiffCriterion']));
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'ModelFitting/Fits/', subjIs{i},'/Fits_locResp_conditioned_unityJdg_jointly']));
    for j = 1:lenMT
        for k = 1:lenDS
            %load data
            C = load(['ModelFitting_updatePrior_', modelTypes{j}{1},'_',...
                modelTypes{j}{2}, strategy{k},'_sub', num2str(subjNs(i)),...
                '.mat'], 'ModelFitting');
            estimatedP{i}{j,k} = C.ModelFitting{end}.P_f;
            minNLL(i,j,k)      = C.ModelFitting{end}.minNLL_f;
            %estimated common-cause prior
            pC1_all(i,j,k,1,:) = estimatedP{i}{j,k}(idx_pC1{j}(1,:)); %cong condition
            pC1_all(i,j,k,2,:) = estimatedP{i}{j,k}(idx_pC1{j}(2,:)); %incong condition
            %estimated criteria
            crt_all(i,j,k,1,:) = estimatedP{i}{j,k}(idx_crt{j}(1,:)); %cong condition
            crt_all(i,j,k,2,:) = estimatedP{i}{j,k}(idx_crt{j}(2,:)); %incong condition
        end
    end
    deltaNLL(i,:,:) = minNLL(i,:,:) - min(min(minNLL(i,:,:)));
end

%% plot AIC
green           = [76,153,0]./255;
colorMap        = [linspace(green(1),1,255)', linspace(green(2),1,255)',...
                    linspace(green(3),1,255)'];
getAIC          = @(nLL, k) 2.*nLL+ 2.*k;
numParams       = [16,17,17,16,17; 14,15,15,14,15; 13,14,14,13,14; 12,13,13,12,13; 11,12,12,11,12]; 
%numTotalTrials  = 320 + 240 + 320*2*2; 
%initialize
[AIC, deltaAIC] = deal(NaN(nS, lenMT, lenDS));
subjI_bestM     = cell(1,nS); %best model name
subjI_bestM_sub = NaN(nS, 2); %subscripts
%the pCommon, sigma_AV_A, sigma_AV_V given the best-fitting decision strategy
[pC1_bestM,pC1_saturated] = deal(NaN(nS, lenC, lenP)); 
[crt_bestM,crt_saturated] = deal(NaN(nS, lenC, lenP)); 
[crt_diff_bestM, crt_diff_saturated] = deal(NaN(nS, lenC));       
bool_plotAIC            = 1;
for i = 1:nS
    %calculate delta AIC for each subject
    AIC(i,:,:)           = getAIC(squeeze(minNLL(i,:,:)), numParams);
    deltaAIC(i,:,:)      = AIC(i,:,:) - min(min(AIC(i,:,:)));
    %find the best model (the model that corresponds to 0 delta AIC)
    idx_bestM            = find(squeeze(deltaAIC(i,:,:)) == 0);
    [row, col]           = ind2sub([lenMT, lenDS], idx_bestM);
    subjI_bestM{i}       = [subjIs{i},modelTypes{row}{1},'-',modelTypes{row}{2},'-', strategy{col}];
    subjI_bestM_sub(i,:) = [row, col];
    disp(subjI_bestM{i}); %display it
    %get the pCommon for all participants
    pC1_bestM(i,:,:)  = squeeze(pC1_all(i,row,col,:,:)); %4 pC1
    %get the sigma for all participants
    crt_bestM(i,1,:)  = estimatedP{i}{row,col}(idx_crt{row}(1,:)); %cong
    crt_bestM(i,2,:)  = estimatedP{i}{row,col}(idx_crt{row}(2,:)); %incong
    crt_diff_bestM(i,:) = squeeze((crt_bestM(i,:,2) - crt_bestM(i,:,1))); 
    
    %assuming 4 crt, 4 pCommon
    [~, idx_temp]         = min(squeeze(deltaAIC(i,1,:)));
    %get the pCommon for all participants
    pC1_saturated(i,:,:)  = squeeze(pC1_all(i,1,idx_temp,:,:)); %4 pC1
    %get the sigma for all participants
    crt_saturated(i,1,:)  = estimatedP{i}{1,idx_temp}(idx_crt{1}(1,:)); %cong
    crt_saturated(i,2,:)  = estimatedP{i}{1,idx_temp}(idx_crt{1}(2,:)); %incong
    crt_diff_saturated(i,:) = squeeze((crt_saturated(i,:,2) - crt_saturated(i,:,1))); 
    %plot
    if bool_plotAIC == 1
        figure
        heatmap(round(squeeze(deltaAIC(i,:,:)),1),'ColorbarVisible', 'on',...
            'YLabel', sprintf('Assumption about \n p_{C=1} and criterion'),...
            'XLabel', sprintf('Assumption about strategies'), 'YData',...
            {'Full - 4 crt', 'Null - 4 crt', 'Full - 1 crt', 'Restricted - 1 crt',...
            'Null - 1 crt'}, 'XData', ...
            {sprintf('MA-posterior_{C=1}'), sprintf('MA-measurements'),...
            sprintf('MA-estimates'),sprintf('MS-posterior_{C=1}'),...
            sprintf('MS-measurements')},'GridVisible','off','Colormap',colorMap); 
        caxis([0, 10]); set(gca,'FontSize',18); 
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.57, 0.4]);
        set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
        %saveas(gcf, ['ModelComparison_deltaAIC_',subjIs{i}], 'pdf'); 
    end
end
pC1_diff_bestM = squeeze(pC1_bestM(:,:,end) - pC1_bestM(:,:,1));
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
    sprintf('Assumption about \n p_{C=1} and criterion'),'XLabel', ...
    sprintf('Assumption about strategies'), 'YData',...
    {'Full - 4 crt', 'Null - 4 crt', 'Full - 1 crt', 'Restricted - 1 crt',...
    'Null - 1 crt'}, 'XData', lbls,'GridVisible','off',...
    'Colormap',colorMap); 
caxis([0, nS]); set(gca,'FontSize',18); %25
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.57, 0.4]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, ['ModelComparison_diffSigma_counts_bestModel_DeltaAICthres',num2str(AIC_thres)], 'pdf'); 

%% plotting the change of sigma and pCommon
plot_bestfitPC1_diffCriterion(pC1_diff_saturated, crt_diff_saturated,...
    deltaAIC, subjI_bestM_sub, 'uni', 1)

%%
plot_bestfitPC1_diffCriterion(pC1_diff_bestM, crt_diff_bestM,...
    deltaAIC, subjI_bestM_sub, 'uni', 1)


