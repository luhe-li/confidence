%This script summarizes the model fits based on the full model, which
%assumes the common-cause prior is different across conditions and phases,
%given the unity judgments and localization responses SEPARATELY. The goal
%was mainly to examine whether the estimated common-cause prior is
%consistent based on the two datasets.
clear all; close all; clc
%subject information
subjNs         = [3,4,5,6,8,9,11,12,13,15,16,17,18,19,20];
lenS           = length(subjNs);
subjIs         = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM',...
                    'SX','ZL','RE','MM'};
%model information
modelTypes     = 'fullModel'; %we only fitted the full model
numDs          = 2; %number of datasets
strat_unityJdg = {'strategyUnity_posteriorC1', ...
                  'strategyUnity_measurements',...
                  'strategyUnity_MAP'}; %3 strategies for unity judgments
len_strat_uj   = length(strat_unityJdg);
strat_locResp  = {'strategyMAP_MA', ...
                  'strategyMAP_MS'}; %2 strategies for localization resps
len_strat_lr   = length(strat_locResp);
cond           = {'cong', 'incong'};
lenC           = length(cond);
phase          = {'pre', 'post'};
lenP           = length(phase);
idx_pC1        = {[8:11;8:11;8:11],[7:10;7:10]}; %indices for the estimated common-cause prior

%% load files and extract the information we want
%initialization
%minimum negative log likelihood
minNLL_unityJdg    = NaN(lenS, len_strat_uj); 
minNLL_locResp     = NaN(lenS, len_strat_lr);
%the common-cause prior estimated based on the unity judgments
pC1_unityJdg       = NaN(lenS, len_strat_uj, lenC, lenP); 
%the difference of the common-cause prior between the pre- and the
%post-learning phase based on the data from the unity judgments
pC1_unityJdg_diff  = NaN(lenS, len_strat_uj, lenC); 
pC1_locResp        = NaN(lenS, len_strat_lr, lenC, lenP);
%the difference of the common-cause prior between the pre- and the
%post-learning phase based on the data from the localization responses
pC1_locResp_diff   = NaN(lenS, len_strat_lr, lenC);
pC1_bestStrat      = NaN(lenS, numDs, lenC, lenP);
%the difference of the common-cause prior based on the best strategy
%2nd dimension: unity judgment, localization responses
pC1_bestStrat_diff = NaN(lenS, numDs, lenC);
%only save the first 6 estimated parameters
numP               = 6;
estP_bestM         = NaN(lenS, numDs, numP);
%p(1)-p(2)  : a_A, b_A, terms that define the relative auditory biases
%p(3)-p(4)  : sigma_A, sigma_V for unimodal trials
%p(5)-p(6)  : sigma_AV_A, sigma_AV_V for bimodal trials

for i = 1:lenS
    %add path
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'ModelFitting/Fits/', subjIs{i},'/Fits_locResp_unityJdg_separately']));
    estP_unityJdg = cell(1,len_strat_uj);
    estP_locResp  = cell(1,len_strat_lr);
    
    %load the model fits given the unity judgments
    for j = 1:len_strat_uj
        C = load(['ModelFitting_updatePrior_', modelTypes, '_', strat_unityJdg{j},...
            '_sub', num2str(subjNs(i)),'.mat'], 'ModelFitting');
        %minimum negative log likelihood
        minNLL_unityJdg(i,j)     = C.ModelFitting{end}.minNLL_f;
        %1st row: pC1_pre_cong,   pC1_post_cong
        %2nd row: pC1_pre_incong, pC1_post_incong
        pC1_temp                 = reshape(C.ModelFitting{end}.P_f(idx_pC1{1}(j,:)),...
                                        [length(phase), length(cond)])';
        estP_unityJdg{j}         = C.ModelFitting{end}.P_f;
        pC1_unityJdg(i,j,:,:)    = pC1_temp;
        pC1_unityJdg_diff(i,j,:) = pC1_temp(:,end) - pC1_temp(:,1);
    end
    %find the index that corresponds to the best strategy
    [~,idx]                   = min(minNLL_unityJdg(i,:));
    pC1_bestStrat(i,1,:,:)    = pC1_unityJdg(i,idx,:,:);
    %disp(strat_unityJdg{idx});
    pC1_bestStrat_diff(i,1,:) = pC1_bestStrat(i,1,:,end) - pC1_bestStrat(i,1,:,1);
    estP_bestM(i,1,:)         = estP_unityJdg{idx}(1:6);
    
    %load the model fits given the localization responses
    for j = 1:len_strat_lr
        C = load(['ModelFitting_updatePrior_', modelTypes, '_', strat_locResp{j},...
            '_sub', num2str(subjNs(i)),'.mat'], 'ModelFitting');
        minNLL_locResp(i,j)     = C.ModelFitting{end}.minNLL_f;
        pC1_temp                = reshape(C.ModelFitting{end}.P_f(idx_pC1{2}(j,:)),...
                                        [length(phase), length(cond)])';
        estP_locResp{j}         = C.ModelFitting{end}.P_f;
        pC1_locResp(i,j,:,:)    = pC1_temp;
        pC1_locResp_diff(i,j,:) = pC1_temp(:,end) - pC1_temp(:,1);
    end
    [~,idx]                   = min(minNLL_locResp(i,:));
    pC1_bestStrat(i,2,:,:)    = pC1_locResp(i,idx,:,:);
    %disp(strat_locResp{idx});
    pC1_bestStrat_diff(i,2,:) = pC1_bestStrat(i,2,:,end) - pC1_bestStrat(i,2,:,1);
    estP_bestM(i,2,:)         = estP_locResp{idx}(1:6);
end

%% plot the estimated common-cause prior based on the unity judgments and 
%the localization responses respectively
bool_save  = 0;
for i = 1:numDs
    pC1_slc     = squeeze(pC1_bestStrat(:,i,:,:)); %selected
    pC1_cmp_slc = squeeze(pC1_bestStrat_diff(:,i,:));
    plot_bestfitPC1(pC1_slc, pC1_slc, pC1_slc, pC1_cmp_slc, pC1_cmp_slc,...
        pC1_cmp_slc, 'rainbow', subjIs, bool_save)
end

%% scatter plot the estimated common-cause prior against the other free parameters
x_lbs = {'a_A','b_A','\sigma_A','\sigma_V','\sigma_{AV,A}','\sigma_{AV,V}'};
figure
for i = 1:numP %6 free parameters
    %the first row is for the congruent condition
    subplot(lenC,numP,i) 
    for j = 1:numDs
        scatter(squeeze(estP_bestM(:,j,i)), squeeze(pC1_bestStrat_diff(:,j,1)),...
            100,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); hold on 
    end
    xlabel(x_lbs{i}); 
    if i == 1
        legend({'Data: unityJdg', 'Data: locResp'}); 
        ylabel('p_{C=1,post} - p_{C=1,pre}');
    end
    set(gca,'FontSize',15);
    
    %the second row is for the incongruent condition
    subplot(lenC,numP,i+numP) 
    for j = 1:numDs
        scatter(squeeze(estP_bestM(:,j,i)), squeeze(pC1_bestStrat_diff(:,j,2)),...
            100,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); hold on 
    end
    xlabel(x_lbs{i}); 
    if i == 1; ylabel('p_{C=1,post} - p_{C=1,pre}'); end
    set(gca,'FontSize',15);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.6]);
set(gcf,'PaperUnits','centimeters','PaperSize',[60 25]);
saveas(gcf, 'freeParamters_diffPC1_scatter', 'pdf')
