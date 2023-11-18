%% load data
%best model:
%'PW':full-MA-MAP,       'SW':samePC1cond-MA-posterior, 'HL':samePC1cond-MS-posterior,
%'YZ':full-MA-MAP,       'NH':samePC1pre-MA-MAP,        'ZZ':full-MA-MAP,
%'BB':full-MS-posterior, 'ZY':full-MA-MAP,              'MR':samePC1pre-MA-MAP,
%'AD':full-MA-MAP,       'SM':full-MA-MAP,              'SX':samePC1pre-MA-m
%'ZL':samePC1pre-MA-MAP  'RE':samePC1cond-MS-posterior  'MM':full-MA-MAP
clear all; close all; clc
subjNs     = [3,4,5,6,8,9,11,12,13,15,16,17,18,19,20];
subjIs     = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL','RE','MM'};
% modelTypes = {'fullModel'};
% strategy   = {'strategyMAP_HA_strategyUnity_MAP'};
modelTypes = {'fullModel','samePC1pre', 'samePC1cond'};
strategy   = {'strategyMAP_MA_strategyUnity_posteriorC1',...
              'strategyMAP_MA_strategyUnity_measurements',...
              'strategyMAP_MA_strategyUnity_MAP',...
              'strategyMAP_MS_strategyUnity_posteriorC1',...
              'strategyMAP_MS_strategyUnity_measurements',...
              'strategyMAP_HA_strategyUnity_MAP'};
cond       = {'cong', 'incong'};
phase      = {'pre', 'post'};
modality   = {'Auditory','Visual'};
bool_data  = [1,1]; %including unity judgments, and localization responses
bool_plot  = [0,0,0]; %1. unity judgment; 2. locResp (raw data); 3. VE
bool_save  = [0,0,0];
idx_mat    = {[8,9;10,11],[8,9;8,10],[8,8;9,9]};

%define global variables and initialize matrices
global sN sI bool_inc lenC lenP lenM lenMT lenDS lenS lenD
lenC  = length(cond); 
lenP  = length(phase); 
lenM  = length(modality);
lenMT = length(modelTypes);
lenDS = length(strategy);
lenS  = 4; %4 stimulus location (A or V)
lenD  = 7; %7 different spatial discrepancy
nS                 = length(subjNs);
[minNLL, deltaNLL] = deal(NaN(nS, lenMT, lenDS));
estimatedP         = cell(1,nS);
pC1_all            = NaN(nS, lenMT, lenDS, lenC, lenP);
for i = 1:nS
    sN = subjNs(i); sI = subjIs{i}; bool_inc = bool_data;
    [estimatedP{i}, minNLL(i,:,:)] = plotModelFits_diffStrategies(modelTypes,...
        strategy, bool_plot, bool_save);
    for j = 1:lenMT
        for k = 1:lenDS
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
getDeviance     = @(nLL_simple, nLL_complex) 2.*(nLL_simple - nLL_complex);
numParams       = [13,14,14,13,14,15; 12,13,13,12,13,14; 11,12,12,11,12,13]; 
%numTotalTrials  = 320 + 240 + 320*2*2; 
%nested hypothesis test
chi_thres_df1   = chi2inv(0.95,1);
chi_thres_df2   = chi2inv(0.95,2);
%initialize
[AIC, deltaAIC] = deal(NaN(nS, lenMT, lenDS));
nestHyp         = NaN(nS, 2); %2 pairs of comparison: null vs full; full vs restricted
subjI_bestM     = cell(1,nS);
subjI_bestM_sub = NaN(nS, 2); %subscripts
%the pCommon, sigma_AV_A, sigma_AV_V given the best-fitting decision strategy
pC1_bestStrategy        = NaN(nS, lenC, lenP); 
weight_shatA            = NaN(1, nS);
bool_plotAIC            = 0;
for i = 1:nS
    %calculate delta AIC for each subject
    AIC(i,:,:)           = getAIC(squeeze(minNLL(i,:,:)), numParams);
    deltaAIC(i,:,:)      = AIC(i,:,:) - min(min(AIC(i,:,:)));
    %find the best model (the model that corresponds to 0 delta AIC)
    idx_bestM            = find(squeeze(deltaAIC(i,:,:)) == 0);
    [row, col]           = ind2sub([lenMT, lenDS], idx_bestM);
    subjI_bestM{i}       = [subjIs{i},'-', modelTypes{row},'-', strategy{col}];
    subjI_bestM_sub(i,:) = [row, col];
    disp(subjI_bestM{i}); %display it

    %get the pCommon for all participants
    pC1_bestStrategy(i,:,:)  = squeeze(pC1_all(i,row,col,:,:));
    %get the sigma for all participants
    weight_shatA(i)      = estimatedP{i}{row, end}(end); 
    
    if bool_plotAIC
        figure
        heatmap(round(squeeze(deltaAIC(i,:,:)),1),'ColorbarVisible', 'on',...
            'YLabel', sprintf('Assumption about p_{C=1}'),...
            'XLabel', sprintf('Assumption about strategies'), 'YData',...
            {'Full', 'Restricted', 'Null'}, 'XData', ...
            {sprintf('MA-posterior_{C=1}'), sprintf('MA-measurements'),...
            sprintf('MA-estimates'),sprintf('MS-posterior_{C=1}'),...
            sprintf('MS-measurements'), sprintf('HA-estimates')},'GridVisible',...
            'off','Colormap',colorMap); 
        caxis([0, 10]); set(gca,'FontSize',18); 
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.35]);
        set(gcf,'PaperUnits','centimeters','PaperSize',[45 12]);
        %saveas(gcf, ['ModelComparison_deltaAIC_',subjIs{i},'_wHA'], 'pdf'); 
    end
    
    %nested hypothesis test
    dev_restricted_vs_full = getDeviance(squeeze(minNLL(i,2,col)), squeeze(minNLL(i,1,col)));
    dev_null_vs_full       = getDeviance(squeeze(minNLL(i,3,col)), squeeze(minNLL(i,1,col)));
    nestHyp(i,1) = (dev_restricted_vs_full > chi_thres_df1);
    nestHyp(i,2) = (dev_null_vs_full       > chi_thres_df2);
end
pC1_diff_bestStrategy = squeeze(pC1_bestStrategy(:,:,end) - pC1_bestStrategy(:,:,1));
%results of nested hypothesis
numSubj_noEffect = sum(nestHyp(:,1)==1 | nestHyp(:,2) ==1);
disp('Nested hypothesis test results: ');
disp(['# of (full > restricted) or (full > null) = ', num2str(numSubj_noEffect)]);

%% count best-fit model
AIC_thres   = 2;
counts_M    = squeeze(sum(deltaAIC < AIC_thres,1));
org         = [255,140,0]./255;
colorMap    = [linspace(1,org(1),255)',linspace(1,org(2),255)',linspace(1,org(3),255)'];
lbls        = {sprintf('MA-posterior_{C=1}'), sprintf('MA-measurements'),...
               sprintf('MA-estimates'),sprintf('MS-posterior_{C=1}'),...
               sprintf('MS-measurements'), sprintf('HA-estimates')};
           
figure
heatmap(counts_M,'ColorbarVisible', 'on','YLabel', ...
    sprintf('Assumption about p_{C=1}'),'XLabel', ...
    sprintf('Assumption about strategies'), 'YData',...
    {'Full', 'Restricted', 'Null'}, 'XData', lbls,'GridVisible','off',...
    'Colormap',colorMap); 
caxis([0, nS]); set(gca,'FontSize',18); %25
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.35]);
set(gcf,'PaperUnits','centimeters','PaperSize',[45 12]);
saveas(gcf, ['ModelComparison_wHA_counts_bestModel_DeltaAICthres',num2str(AIC_thres)], 'pdf'); 

%% plotting the change of sigma and pCommon
plot_bestfitPC1_diffSigma(pC1_diff_bestStrategy, sigma_diff_bestStrategy,...
    deltaAIC, subjI_bestM_sub, 'uni', 1)


