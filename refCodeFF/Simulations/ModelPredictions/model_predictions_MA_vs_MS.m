%% This script shows the difference of model predictions by 'full',
%'restricted', and 'null' models using simulated data
clear all; close all; clc; rng(2);
addpath(genpath(['/Users/Guinevere/Desktop/NYU/Project1/Experiment code_project1/',...
                        'ModelFitting/bads-master']));
strat_locResp        = {'MA','MS'};
model.numSD          = 5;
model.modality       = {'A','V'};
model.strategy_unity = 'posteriorC1';
model.numBins_A      = 100;
model.numBins_V      = 100;
model.bins_r         = -35:0.1:30;
model.numBins_r      = length(model.bins_r);

%define parameters
pC1                  = 0.5;
s_A                  = -12:4:12;
s_V                  = -12:4:12;
AV_allcombs          = combvec(s_A,s_V);
spatialD             = unique(AV_allcombs(2,:) - AV_allcombs(1,:));
a_A                  = 1;
b_A                  = 0; 
sigma_AV_A           = 10;
sigma_AV_V           = 2.5;
epsilon              = NaN; 
mu_P                 = 0;
sigma_P              = 30;
lapse                = 0.03;
sigma_r              = 1;

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
global lenD lenM lenS
lenM = length(model.modality);
lenD = length(spatialD); 
lenS = length(s_A);
%initialize a matrix
probResp_true = NaN(length(strat_locResp),lenS, lenS, lenM, model.numBins_r);
for m = 1:length(strat_locResp)
    model.strategy_MAP = strat_locResp{m};
    disp(m)
    %get model prediction
    [~,R] = get_mPred(pC1,epsilon,s_A, s_V, a_A, b_A,...
        sigma_AV_A,sigma_AV_V, CI, mu_P, lapse, model, []); 
    %the last argument is for data (if we are not calculating nLL, leave it [])
    probResp_true(m,:,:,:,:,:) = predict_prob_r(R.MAP, ...
        R.p_mAmV_given_sAsV,sigma_r, model);
end

%% simulate data given ground truth 
%visualize 
selected_A_loc  = 1;
for m = 1:length(strat_locResp)
    probResp_true_ijk = squeeze(probResp_true(m,selected_A_loc,:,:,:));
    plt_locResp_sim(probResp_true_ijk, model.bins_r, model, s_V, strat_locResp{m},1);
end



