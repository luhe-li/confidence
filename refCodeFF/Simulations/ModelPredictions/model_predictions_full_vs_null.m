%% This script shows the difference of model predictions by 'full',
%'restricted', and 'null' models using simulated data
clear all; close all; clc; rng(2);
% addpath(genpath(['/Users/Guinevere/Desktop/NYU/Project1/Experiment code_project1/',...
%                         'ModelFitting/bads-master']));
modelType            = {'full','restricted','null'};
model.cond           = {'cong', 'incong'};
model.phase          = {'pre', 'post'};
model.numSD          = 5;
model.modality       = {'A','V'};
model.strategy_MAP   = 'MA';
model.strategy_unity = 'MAP';
model.numBins_A      = 100;
model.numBins_V      = 100;
model.bins_r         = -30:0.1:30;
model.numBins_r      = length(model.bins_r);

%define parameters
pC1                  = {[0.2,0.9;0.6,0.1],[0.5,0.9;0.5,0.1],[0.6,0.6;0.2,0.2]};
s_A                  = -12:4:12;
s_V                  = -12:4:12;
AV_allcombs          = combvec(s_A,s_V);
spatialD             = unique(AV_allcombs(2,:) - AV_allcombs(1,:));
a_A                  = 1;
b_A                  = 0; 
sigma_AV_A           = 5;
sigma_AV_V           = 2.5;
epsilon              = 3; 
mu_P                 = 0;
sigma_P              = 100;
lapse                = 0.03;
sigma_r              = 1.5;
nT_perPair           = 20; %how many unity judgment per AV pair

%constants for fitting the unity judgment
CI.J_A               = sigma_AV_A^2;
CI.J_V               = sigma_AV_V^2;
CI.J_P               = sigma_P^2;
CI.constC1           = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
CI.constC1_shat      = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
CI.constC2_1         = CI.J_A + CI.J_P;
CI.constC2_1_shat    = 1/CI.J_A + 1/CI.J_P;
CI.constC2_2         = CI.J_V + CI.J_P;
CI.constC2_2_shat    = 1/CI.J_V + 1/CI.J_P; 

%% Given the true pCommon for both conditions, we first get the model 
%predictions on the proportion of reporting 'common cause' for each pair of
%audiovisual stimulus

%define global variables
global lenC lenD lenP lenM lenS
lenC = length(model.cond); 
lenP = length(model.phase); 
lenM = length(model.modality);
lenD = length(spatialD); 
lenS = length(s_A);
%initialize a matrix
[uniJdg_perc_true, uniJdg_perc_sim] = deal(NaN(length(modelType), lenC, lenP, lenS, lenS));
uniJdg_binary_sim = NaN(length(modelType), lenC, lenP, lenS, lenS, nT_perPair);
[pC1_resp_lenD_true, uniJdg_perc_sim_lenD,uniJdg_perc_sim_SD_lenD] = ...
    deal(NaN(length(modelType), lenC, lenP, lenD));
for m = 1:length(modelType)
    for i = 1:lenC
        for j = 1:lenP
            %get model prediction
            [~,R] = get_mPred(pC1{m}(i,j),epsilon,s_A, s_V, a_A, b_A, sigma_AV_A,...
                sigma_AV_V, CI, mu_P, lapse, model, []); 
            %the last argument is for data (if we are not calculating nLL, leave it [])
            uniJdg_perc_true(m,i,j,:,:) = R.pC1_given_sAsV;
        end
    end
    % simulate data given ground truth 
    uniJdg_perc_sim(m,:,:,:,:) = sim_reportingC1(squeeze(uniJdg_perc_true(m,:,:,:,:)), nT_perPair);
    %simulate binary responses given a p(reporting C=1)
    %e.g., if p=0.5 and nT_perPair = 4, then sim_binaryJdg gives [1,1,0,0]
    uniJdg_binary_sim(m,:,:,:,:,:) = sim_binaryJdg(squeeze(uniJdg_perc_sim(m,:,:,:,:)), nT_perPair);

    %convert the data from (lenS x lenS) to (lenS x lenD)
    [~, pC1_resp_lenD_true(m,:,:,:), ~, numT_AV] = ...
        reshapeResp(squeeze(uniJdg_perc_true(m,:,:,:,:)),nT_perPair);
    [~, uniJdg_perc_sim_lenD(m,:,:,:), uniJdg_perc_sim_SD_lenD(m,:,:,:), ~] = ...
        reshapeResp(squeeze(uniJdg_perc_sim(m,:,:,:,:)),nT_perPair);
end

%% plot
for m = 1:length(modelType)
%     plt_uniJdg_perc_sim(spatialD, squeeze(pC1_resp_lenD_true(m,:,:,:)),[],...
%         squeeze(uniJdg_perc_sim_lenD(m,:,:,:)),squeeze(uniJdg_perc_sim_SD_lenD(m,:,:,:)),...
%         numT_AV, modelType{m},0);
    plt_uniJdg_perc_sim(spatialD, squeeze(pC1_resp_lenD_true(m,:,:,:)),[],...
        [],[],numT_AV, modelType{m},1);
end

%% fit models (with different decision strategies) to the simulated data
%%first initialize matrices)
% fit_pCommon_strat         = cell(1,length(strat_unity));
% pC1_resp_mPred_strat      = NaN(length(strat_unity), lenC, lenP, lenS, lenS);
% pC1_resp_lenD_mPred_strat = NaN(length(strat_unity), lenC, lenP, lenD); 
% 
% for s = 1:length(strat_unity)
%     model.strategy_unity = strat_unity{s};
%     %fit the model to the similated data
%     fit_pCommon_strat{s} = fit_uniJdg_sim_strat(uniJdg_binary_sim, s_A, s_V,...
%         a_A, b_A, sigma_AV_A, sigma_AV_V, CI, mu_P, lapse, model);
% 
%     %get model predictions given the best-fit pCommon
%     fit_pCommon_temp = fit_pCommon_strat{s}(1:lenC*lenP);
%     epsilon_temp = fit_pCommon_strat{s}(end);
%     for i = 1:lenC
%         for j = 1:lenP
%             [~,R] = get_mPred(fit_pCommon_temp((i-1)*lenC+j), epsilon_temp,...
%                 s_A, s_V, a_A, b_A, sigma_AV_A, sigma_AV_V, CI, mu_P, lapse,...
%             model, []);
%             pC1_resp_mPred_strat(s,i,j,:,:) = R.pC1_given_sAsV;
%         end
%     end
%     %reshape
%     [~, pC1_resp_lenD_mPred_strat(s,:,:,:), ~, ~] = reshapeResp(squeeze(...
%         pC1_resp_mPred_strat(s,:,:,:,:)), nT_perPair);
% end
% 
% for s = 1:length(strat_unity)
%     %plot the simulated data with model 
%     %pC1_resp = pC1_resp_lenD_true;
%     if s == 1; pC1_resp = pC1_resp_lenD_true; else; pC1_resp = [];end
%     plt_uniJdg_perc_sim(spatialD, pC1_resp, squeeze(pC1_resp_lenD_mPred_strat(s,:,:,:)),...
%         uniJdg_perc_sim_lenD, uniJdg_perc_sim_SD_lenD, numT_AV, strat_unity{s},0);
% end






